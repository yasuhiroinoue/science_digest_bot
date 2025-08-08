import os
import re
import time
import json
import hashlib
import quopri
import requests
import feedparser
import xml.etree.ElementTree as ET
from datetime import datetime
from dotenv import load_dotenv
from dateutil.parser import parse as parse_date
from google import genai
from google.genai import types

"""
Science RSS -> Abstract fetch (Crossref/PubMed) -> Gemini要約(日本語) -> Discord Webhook 投稿
- Quoted-Printable(QP)を自動デコード
- 文字列の再パース(「Key: Value」)を廃止し、辞書から直接Markdownを生成
- 再試行つきWebhook送信
"""

# =====================
# 環境変数の読み込み
# =====================
load_dotenv()

region = os.getenv("REGION", "global")
project_id = os.getenv("PROJECT_ID", "default_project_id")
rss_url = os.getenv("RSS_URL", "default_rss_url")
webhook_urls = [u.strip() for u in os.getenv("WEBHOOK_URLS", "default_webhook_urls").split(",") if u.strip()]
hash_file_path = os.getenv("HASH_FILE_PATH", "latest_entry_hash.txt")

HEADERS = {
    "User-Agent": os.getenv(
        "USER_AGENT",
        "DOIScraperBot/1.0 (mailto:your-email@example.com)",
    )
}

MODEL_ID = os.getenv("MODEL_ID", "gemini-2.5-flash-lite-preview-06-17")
try:
    BUDGET = float(os.getenv("BUDGET", "0"))
except ValueError:
    BUDGET = 0.0

# =====================
# Gemini 設定
# =====================

generate_content_config = types.GenerateContentConfig(
    temperature=1,
    top_p=0.95,
    max_output_tokens=8192,
    system_instruction=
    """
あなたは高度な英語-日本語翻訳者です。提供された英語の論文要旨を学術的な聴衆に向けて日本語に要約してください。要約は100文字以内に抑え、論文の核心を簡潔に表現してください。文字数の制限は論文の内容の複雑さによって柔軟に調整可能です。要約が不可能な場合、または「Abstract not found.」と入力された場合は、「要約なし」と返してください。論文要旨が短い場合、直接的な日本語訳を提供してください。プロセスの結果として、要約または翻訳されたテキストのみを返すようにしてください。「要旨:」という見出しを付けることは禁止します。
""",
    thinking_config=types.ThinkingConfig(thinking_budget=BUDGET),
    response_mime_type="text/plain",
)

gemini_client = genai.Client(vertexai=True, project=project_id, location=region)


def generate_text(prompt: str) -> str:
    response = gemini_client.models.generate_content(
        model=MODEL_ID, contents=prompt, config=generate_content_config
    )
    return response.text


# =====================
# RSS ハッシュ管理
# =====================

def get_latest_entry_hash(rss_data):
    latest_entry = rss_data["entries"][0]
    unique_str = latest_entry.get("title", "") + latest_entry.get("published", "")
    return hashlib.md5(unique_str.encode()).hexdigest()


def save_latest_entry_hash(hash_str, path):
    with open(path, "w") as f:
        f.write(hash_str)


def read_latest_entry_hash(path):
    if os.path.exists(path):
        with open(path, "r") as f:
            return f.read()
    return None


# =====================
# RSS 取得
# =====================

def get_latest_articles_from_rss(feed_url, hash_file_path):
    feed = feedparser.parse(feed_url)
    current_hash = get_latest_entry_hash(feed)
    previous_hash = read_latest_entry_hash(hash_file_path)

    if current_hash == previous_hash:
        print("No new entries, skipping processing.")
        return [], False

    articles = []
    for entry in feed.entries:
        doi = entry.get("dc_identifier", "").split("doi:")[-1].strip()
        title = entry.title
        updated = entry.get("updated", "No date available")
        parsed_date = parse_date(updated) if updated != "No date available" else "Date unavailable"
        authors = entry.get("author", "Authors not listed")
        url = entry.link
        articles.append({
            "DOI": doi,
            "Title": title,
            "Authors": authors,
            "Date": parsed_date,
            "URL": url,
        })
    return articles, True


# =====================
# Abstract 取得 (Crossref / PubMed)
# =====================

def get_crossref_abstract(doi: str, headers: dict) -> str | None:
    """Crossref からアブストラクト取得。HTMLタグはここでは除去しない。"""
    api_url = f"https://api.crossref.org/works/{doi}"
    print(f"Requesting Crossref: {api_url}")
    try:
        response = requests.get(api_url, headers=headers, timeout=10)
        response.raise_for_status()
        try:
            data = response.json()
        except json.JSONDecodeError:
            print(f"Error: Failed to decode JSON response from Crossref for DOI {doi}.")
            return None
        message = data.get("message")
        if not message:
            print("Could not find 'message' key in Crossref response.")
            return None
        abstract = message.get("abstract")
        if abstract:
            print(f"Abstract found via Crossref for DOI {doi}.")
            return abstract
        print(f"Abstract not found in Crossref metadata for DOI {doi}.")
        return None
    except requests.exceptions.Timeout:
        print(f"Error: Request timed out fetching Crossref for DOI {doi}.")
        return None
    except requests.exceptions.HTTPError as http_err:
        print(f"Crossref HTTP error for DOI {doi}: {http_err} - Status: {response.status_code}")
        return None
    except requests.exceptions.RequestException as req_err:
        print(f"Crossref request error for DOI {doi}: {req_err}")
        return None
    except Exception as e:
        print(f"Unexpected error during Crossref fetch for DOI {doi}: {e}")
        return None


def get_pubmed_abstract(doi: str, headers: dict) -> str | None:
    """PubMed E-utilities で DOI -> PMID -> Abstract を取得。"""
    print(f"Attempting PubMed search for DOI: {doi}")

    # --- ESearch: DOI -> PMID ---
    esearch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    email_match = re.search(r"mailto:([^)]+)", headers.get("User-Agent", ""))
    email = email_match.group(1) if email_match else "unknown@example.com"
    esearch_params = {"db": "pubmed", "term": doi, "retmode": "json", "tool": "DOIScraperBot", "email": email}

    pmid = None
    try:
        time.sleep(0.4)  # NCBI rate limits
        response_search = requests.get(esearch_base, params=esearch_params, headers=headers, timeout=15)
        response_search.raise_for_status()
        try:
            search_data = response_search.json()
        except json.JSONDecodeError:
            print(f"Error: Failed to decode JSON from PubMed ESearch for DOI {doi}.")
            return None
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            print(f"No PubMed ID (PMID) found for DOI: {doi}")
            return None
        pmid = id_list[0]
        print(f"Found PMID: {pmid} for DOI: {doi}")
    except requests.exceptions.Timeout:
        print(f"Error: Timeout during PubMed ESearch for DOI {doi}.")
        return None
    except requests.exceptions.HTTPError as http_err:
        print(f"PubMed ESearch HTTP error for DOI {doi}: {http_err} - Status: {response_search.status_code}")
        return None
    except requests.exceptions.RequestException as req_err:
        print(f"PubMed ESearch request error for DOI {doi}: {req_err}")
        return None
    except Exception as e:
        print(f"Unexpected error during PubMed ESearch for DOI {doi}: {e}")
        return None

    # --- EFetch: PMID -> Abstract(XML) ---
    print(f"Fetching abstract from PubMed using PMID: {pmid}")
    efetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    efetch_params = {
        "db": "pubmed",
        "id": pmid,
        "retmode": "xml",
        "rettype": "abstract",
        "tool": "DOIScraperBot",
        "email": email,
    }

    try:
        time.sleep(0.4)
        response_fetch = requests.get(efetch_base, params=efetch_params, headers=headers, timeout=15)
        response_fetch.raise_for_status()
        try:
            xml_root = ET.fromstring(response_fetch.text)
            abstract_elements = xml_root.findall(".//AbstractText")
            if not abstract_elements:
                print(f"Abstract text not found in PubMed XML for PMID: {pmid}")
                return None
            abstract_parts = []
            for elem in abstract_elements:
                text = elem.text
                label = elem.get("Label")
                if text:
                    part = f"{label}: {text.strip()}" if label else text.strip()
                    abstract_parts.append(part)
            full_abstract = "\n".join(abstract_parts)
            if full_abstract:
                print(f"Abstract found via PubMed for PMID: {pmid}.")
                return full_abstract
            print(f"Abstract text empty in PubMed XML for PMID: {pmid}")
            return None
        except ET.ParseError:
            print(f"Error: Failed to parse XML from PubMed EFetch for PMID {pmid}.")
            return None
        except Exception as parse_err:
            print(f"Error during PubMed XML parsing for PMID {pmid}: {parse_err}")
            return None
    except requests.exceptions.Timeout:
        print(f"Error: Timeout during PubMed EFetch for PMID {pmid}.")
        return None
    except requests.exceptions.HTTPError as http_err:
        print(f"PubMed EFetch HTTP error for PMID {pmid}: {http_err} - Status: {response_fetch.status_code}")
        return None
    except requests.exceptions.RequestException as req_err:
        print(f"PubMed EFetch request error for PMID {pmid}: {req_err}")
        return None
    except Exception as e:
        print(f"Unexpected error during PubMed EFetch for PMID {pmid}: {e}")
        return None


# =====================
# テキストユーティリティ
# =====================

def decode_qp_if_needed(s: str) -> str:
    """Quoted-Printableの痕跡があればUTF-8としてデコード"""
    if not isinstance(s, str):
        return s
    if re.search(r"=\r?\n|=[0-9A-F]{2}", s, re.IGNORECASE):
        try:
            return quopri.decodestring(s).decode("utf-8", errors="replace")
        except Exception:
            return s
    return s


def sanitize_text(s: str) -> str:
    return decode_qp_if_needed(s).strip() if isinstance(s, str) else s


def format_date_for_md(d) -> str:
    if isinstance(d, datetime):
        return d.isoformat()
    return str(d)


def format_authors(info):
    authors = str(info.get("Authors", "Authors not listed")).replace("\n", "")
    author_names = [a.strip() for a in authors.split(",") if a.strip()]
    if len(author_names) > 4:
        formatted = f"{author_names[0]}, {author_names[1]}, {author_names[2]} et al."
    else:
        formatted = ", ".join(author_names) if author_names else "Authors not listed"
    return [formatted]


def build_discord_markdown(info: dict, abstract_ja: str) -> str:
    title = sanitize_text(info.get("Title", "(No title)"))
    authors_line = format_authors(info)[0]
    date_str = format_date_for_md(info.get("Date", "Date unavailable"))
    url = info.get("URL", "")
    abstract_block = sanitize_text(abstract_ja) or "要約なし"
    abstract_block = "\n".join([ln for ln in abstract_block.splitlines() if ln.strip()])
    md = (
        f"- **{title}**\n"
        f" {authors_line}\n"
        f" {date_str}\n"
        f" ```{abstract_block}```\n"
        f" {url}\n"
    )
    return md


# =====================
# Discord 投稿
# =====================

def generate_and_send_messages(extracted_info, webhook_urls):
    count = 1
    for info in extracted_info:
        print(count)
        count += 1

        # 英文アブストラクト(または"Abstract not found.")
        abstract_en = str(info.get("summary", "Abstract not found."))
        abstract_en = sanitize_text(abstract_en)

        # Gemini で和訳・要約
        abstract_ja = sanitize_text(generate_text(abstract_en))
        print(abstract_ja)  # debug

        discord_message = build_discord_markdown(info, abstract_ja)
        print(discord_message)  # debug

        payload = {"content": discord_message}
        max_retries = 5
        retry_delay = 1

        for url in webhook_urls:
            for attempt in range(max_retries):
                try:
                    response = requests.post(url, json=payload, timeout=15)
                except Exception as post_err:
                    print(f"Post error: {post_err}")
                    response = None
                if response and response.status_code in (200, 204):
                    break
                else:
                    status = response.status_code if response else "NO_RESP"
                    reason = response.reason if response else "Exception"
                    print(f"Error: {status} - {reason}. Retrying...")
                    time.sleep(retry_delay)
            else:
                print(f"Failed to post after {max_retries} attempts.")

        time.sleep(1)  # rate-limit


# =====================
# メイン
# =====================

def main():
    articles, has_new_entries = get_latest_articles_from_rss(rss_url, hash_file_path)
    if not has_new_entries:
        print("No new RSS entries found.")
        return

    extracted_info = []
    processed_count = 0

    print(f"Found {len(articles)} new articles. Processing all...")
    for article in articles:
        doi = article.get("DOI")
        if not doi:
            print(f"Skipping article with missing DOI: {article.get('Title')}")
            continue

        print(f"\nProcessing article: {article.get('Title')} (DOI: {doi})")
        abstract = None

        # 1) Crossref
        abstract = get_crossref_abstract(doi, HEADERS)

        # 2) Crossrefで無ければPubMed
        if abstract is None:
            print("Abstract not found via CrossRef, trying PubMed...")
            time.sleep(0.5)
            abstract = get_pubmed_abstract(doi, HEADERS)

        # 3) Crossref由来のタグっぽい場合は除去
        if abstract and isinstance(abstract, str):
            if abstract.strip().startswith("<"):
                print("Cleaning potential tags from abstract...")
                cleaned = re.sub("<[^>]*>", "", abstract).strip()
                abstract = cleaned if cleaned else "Abstract found but empty after cleaning."
        if abstract == "Abstract not found.":
            abstract = None

        # 4) なお取れない場合のデフォルト
        if abstract is None:
            print(f"Abstract could not be retrieved for DOI {doi} from any source.")
            article["summary"] = "Abstract not found."
        else:
            print("Successfully obtained abstract.")
            article["summary"] = abstract

        extracted_info.append(article)
        processed_count += 1
        time.sleep(1)

    print(f"\nFinished processing {processed_count} articles.")

    if extracted_info:
        print("Sending messages to Discord...")
        generate_and_send_messages(extracted_info, webhook_urls)
        # 成功後にハッシュを保存
        try:
            save_latest_entry_hash(get_latest_entry_hash(feedparser.parse(rss_url)), hash_file_path)
        except Exception as e:
            print(f"Warning: failed to save latest hash: {e}")


if __name__ == "__main__":
    main()
