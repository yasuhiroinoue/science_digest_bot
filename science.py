import feedparser
import requests
import feedparser
import requests
import hashlib
import os
import time
import json  # Needed for JSON decoding errors
import xml.etree.ElementTree as ET  # Needed for PubMed XML parsing
import re  # Needed for cleaning abstract tags
from dotenv import load_dotenv
from dateutil.parser import parse as parse_date  # To parse date
from google import genai
from google.genai import types


# Load environment variables from a .env file
load_dotenv()

# Replace hard-coded values with values from environment variables
region = os.getenv("REGION", "global")
project_id = os.getenv("PROJECT_ID", "default_project_id")
rss_url = os.getenv("RSS_URL", "default_rss_url")
webhook_urls = os.getenv("WEBHOOK_URLS", "default_webhook_urls").split(",")
hash_file_path = os.getenv("HASH_FILE_PATH", "latest_entry_hash.txt")

# --- Headers for API requests ---
HEADERS = {
    "User-Agent": "DOIScraperBot/1.0 (mailto:your-email@example.com)"
}  # Use email from .env?

# Gemini
MODEL_ID = os.getenv("MODEL_ID", "gemini-2.5-flash-lite-preview-06-17")
BUDGET = os.getenv("BUDGET", 0)

generate_content_config = types.GenerateContentConfig(
    temperature=1,
    top_p=0.95,
    # top_k= 64,
    max_output_tokens=8192,
    system_instruction="""
あなたは高度な英語-日本語翻訳者です。提供された英語の論文要旨を学術的な聴衆に向けて日本語に要約してください。要約は100文字以内に抑え、論文の核心を簡潔に表現してください。文字数の制限は論文の内容の複雑さによって柔軟に調整可能です。要約が不可能な場合、または「Abstract not found.」と入力された場合は、「要約なし」と返してください。論文要旨が短い場合、直接的な日本語訳を提供してください。プロセスの結果として、要約または翻訳されたテキストのみを返すようにしてください。「要旨:」という見出しを付けることは禁止します。
""",
    thinking_config=types.ThinkingConfig(
        thinking_budget=BUDGET,
    ),
    response_mime_type="text/plain",
)

gemini_client = genai.Client(vertexai=True, project=project_id, location=region)


# Gemini response
def generate_text(prompt):
    response = gemini_client.models.generate_content(
        model=MODEL_ID, contents=prompt, config=generate_content_config
    )
    return response.text


def get_latest_entry_hash(rss_data):
    latest_entry = rss_data["entries"][0]
    unique_str = latest_entry.get("title", "") + latest_entry.get("published", "")
    return hashlib.md5(unique_str.encode()).hexdigest()


def save_latest_entry_hash(hash_str, hash_file_path):
    with open(hash_file_path, "w") as file:
        file.write(hash_str)


def read_latest_entry_hash(hash_file_path):
    if os.path.exists(hash_file_path):
        with open(hash_file_path, "r") as file:
            return file.read()
    return None


def get_latest_articles_from_rss(feed_url, hash_file_path):
    feed = feedparser.parse(feed_url)
    current_hash = get_latest_entry_hash(feed)
    previous_hash = read_latest_entry_hash(hash_file_path)

    if current_hash == previous_hash:
        print("No new entries, skipping processing.")
        return [], False  # No articles

    articles = []
    for entry in feed.entries:
        doi = entry.get("dc_identifier", "").split("doi:")[-1].strip()
        title = entry.title
        updated = entry.get(
            "updated", "No date available"
        )  # Get date from the updated key
        parsed_date = (
            parse_date(updated)
            if updated != "No date available"
            else "Date unavailable"
        )  # parse date
        authors = entry.get("author", "Authors not listed")
        url = entry.link
        articles.append(
            {
                "DOI": doi,
                "Title": title,
                "Authors": authors,
                "Date": parsed_date,  # 記事の辞書に日付を追加
                "URL": url,
            }
        )
    return articles, True


# --- Abstract Fetching Functions (CrossRef and PubMed Fallback) ---


def get_crossref_abstract(doi: str, headers: dict) -> str | None:
    """
    Crossref APIを使用して指定されたDOIの論文のアブストラクトを取得します。
    Args:
        doi: 取得したい論文のDOI
        headers: リクエストヘッダー
    Returns:
        アブストラクトの文字列。見つからない場合やエラーの場合はNone。
    """
    api_url = f"https://api.crossref.org/works/{doi}"
    print(f"Requesting Crossref: {api_url}")

    try:
        # Use a reasonable timeout
        response = requests.get(api_url, headers=headers, timeout=10)
        response.raise_for_status()  # Check for HTTP errors

        try:
            data = response.json()
        except json.JSONDecodeError:
            print(f"Error: Failed to decode JSON response from Crossref for DOI {doi}.")
            return None

        message = data.get("message")
        if message:
            abstract = message.get("abstract")
            if abstract:
                print(f"Abstract found via Crossref for DOI {doi}.")
                # Return raw abstract which might contain tags
                return abstract
            else:
                # Abstract key exists but is null/empty, or abstract key doesn't exist
                print(f"Abstract not found in Crossref metadata for DOI {doi}.")
                return None
        else:
            print(
                f"Could not find 'message' key in Crossref API response for DOI {doi}."
            )
            return None

    except requests.exceptions.Timeout:
        print(f"Error: Request timed out fetching Crossref for DOI {doi}.")
        return None
    except requests.exceptions.HTTPError as http_err:
        # Log specific HTTP errors
        print(
            f"Crossref HTTP error for DOI {doi}: {http_err} - Status: {response.status_code}"
        )
        # Don't return error message here, just None, let main handle default
        return None
    except requests.exceptions.RequestException as req_err:
        print(f"Crossref request error for DOI {doi}: {req_err}")
        return None
    except Exception as e:
        # Catch unexpected errors during Crossref fetch
        print(f"Unexpected error during Crossref fetch for DOI {doi}: {e}")
        return None


def get_pubmed_abstract(doi: str, headers: dict) -> str | None:
    """
    PubMed API (E-utilities) を使用してDOIからアブストラクトを取得します。
    Args:
        doi: 取得したい論文のDOI
        headers: リクエストヘッダー
    Returns:
        アブストラクトの文字列。見つからない場合やエラーの場合はNone。
    """
    print(f"Attempting PubMed search for DOI: {doi}")

    # --- 1. ESearch: DOI to PMID ---
    esearch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    # Extract email for NCBI tool/email params if possible
    email_match = re.search(r"mailto:([^)]+)", headers.get("User-Agent", ""))
    email = email_match.group(1) if email_match else "unknown@example.com"
    esearch_params = {
        "db": "pubmed",
        "term": doi,
        "retmode": "json",
        "tool": "DOIScraperBot",
        "email": email,
    }

    pmid = None
    try:
        time.sleep(0.4)  # Adhere to NCBI rate limits (3 requests/sec without API key)
        response_search = requests.get(
            esearch_base, params=esearch_params, headers=headers, timeout=15
        )
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
        print(
            f"PubMed ESearch HTTP error for DOI {doi}: {http_err} - Status: {response_search.status_code}"
        )
        return None
    except requests.exceptions.RequestException as req_err:
        print(f"PubMed ESearch request error for DOI {doi}: {req_err}")
        return None
    except Exception as e:
        print(f"Unexpected error during PubMed ESearch for DOI {doi}: {e}")
        return None

    if not pmid:
        return None  # Exit if no PMID found

    # --- 2. EFetch: PMID to Abstract ---
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
        time.sleep(0.4)  # Adhere to NCBI rate limits
        response_fetch = requests.get(
            efetch_base, params=efetch_params, headers=headers, timeout=15
        )
        response_fetch.raise_for_status()

        # --- 3. XML Parsing ---
        try:
            xml_root = ET.fromstring(response_fetch.text)
            abstract_elements = xml_root.findall(".//AbstractText")
            if not abstract_elements:
                print(f"Abstract text not found in PubMed XML for PMID: {pmid}")
                return None

            # Join structured abstracts (e.g., BACKGROUND:, METHODS:, RESULTS:)
            abstract_parts = []
            for elem in abstract_elements:
                text = elem.text
                label = elem.get("Label")
                if text:
                    # Prepend label if it exists
                    part = f"{label}: {text.strip()}" if label else text.strip()
                    abstract_parts.append(part)

            full_abstract = "\n".join(abstract_parts)
            if full_abstract:
                print(f"Abstract found via PubMed for PMID: {pmid}.")
                # Assume PubMed abstract is clean text
                return full_abstract
            else:
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
        print(
            f"PubMed EFetch HTTP error for PMID {pmid}: {http_err} - Status: {response_fetch.status_code}"
        )
        return None
    except requests.exceptions.RequestException as req_err:
        print(f"PubMed EFetch request error for PMID {pmid}: {req_err}")
        return None
    except Exception as e:
        print(f"Unexpected error during PubMed EFetch for PMID {pmid}: {e}")
        return None


# --- End of Abstract Fetching Functions ---


def format_authors(info):
    authors = info["Authors"].replace("\n", "")  # Remove any newline characters
    author_list = []

    # Split the author names into a list
    author_names = authors.split(", ")

    if len(author_names) > 4:
        # Select first three authors and append 'et al.'
        formatted_authors = (
            f"{author_names[0]}, {author_names[1]}, {author_names[2]} et al."
        )
    else:
        # Join all authors because there are three or less
        formatted_authors = ", ".join(author_names)

    # Append the formatted string to author_list
    author_list.append(formatted_authors)

    return author_list


def format_to_markdown(content_message):
    # Extract the information from the content message
    info = {}
    for line in content_message.split("\n"):
        if line.strip():  # Check if line is not empty
            key, value = line.split(":", 1)
            info[key.strip()] = value.strip()

    # Remove blank lines from the Summary
    if "Summary" in info:
        lines = info["Summary"].split("\n")
        filtered_summary = "\n".join([line for line in lines if line.strip()])
        info["Summary"] = filtered_summary

    # Format the extracted information into markdown
    markdown_format = (
        f"- **{info['Title']}**\n"
        f" {info['Authors']}\n"
        f" {info['Date']}\n"
        f" ```{filtered_summary}```\n"
        f" {info['URL']}\n"
    )

    return markdown_format


def generate_and_send_messages(extracted_info, webhook_urls):
    # client = AnthropicVertex(region=region, project_id=project_id)
    count = 1
    for info in extracted_info:
        print(count)
        count += 1
        content_message = f"Summary: {info['summary']}\n"

        abstract = generate_text(content_message)
        # debug
        print(abstract)

        authors = format_authors(info)
        content = f"Title: {info['Title']}\nAuthors: {authors[0]}\nDate: {info['Date']}\nSummary: {abstract}\nURL: {info['URL']}"

        print(content)

        discord_message = format_to_markdown(content)

        # Post the generated content to Discord
        payload = {"content": discord_message}
        max_retries = 5  # 最大再試行回数
        retry_delay = 1  # 再試行の間隔（秒）

        for url in webhook_urls:
            for attempt in range(max_retries):
                response = requests.post(url, json=payload)

                if response.status_code in [200, 204]:
                    break  # 成功したらループを抜ける
                else:
                    print(
                        f"Error: {response.status_code} - {response.reason}. Retrying..."
                    )
                    print(discord_message)
                    time.sleep(retry_delay)  # 指定された時間待つ
            else:
                # すべての再試行が失敗した場合
                print(f"Failed to post after {max_retries} attempts.")

        time.sleep(1)  # Adjust based on your rate limit policies


def main():
    articles, has_new_entries = get_latest_articles_from_rss(rss_url, hash_file_path)
    if not has_new_entries:
        print("No new RSS entries found.")
        return

    extracted_info = []
    processed_count = 0
    # fetch_limit removed as requested

    print(f"Found {len(articles)} new articles. Processing all...")

    for article in articles:  # Process all articles
        doi = article.get("DOI")
        if not doi:
            print(f"Skipping article with missing DOI: {article.get('Title')}")
            continue

        print(f"\nProcessing article: {article.get('Title')} (DOI: {doi})")
        abstract = None
        source = None

        # 1. Try CrossRef
        abstract = get_crossref_abstract(doi, HEADERS)
        source = "CrossRef"

        # 2. If CrossRef fails, try PubMed
        if abstract is None:
            print("Abstract not found via CrossRef, trying PubMed...")
            time.sleep(0.5)  # Small delay between API calls
            abstract = get_pubmed_abstract(doi, HEADERS)
            source = "PubMed"

        # 3. Clean the abstract if found and potentially contains tags (mainly from CrossRef)
        if abstract and isinstance(abstract, str):
            # Basic cleaning for XML/HTML tags
            if abstract.strip().startswith("<"):
                print("Cleaning potential tags from abstract...")
                cleaned_abstract = re.sub("<[^>]*>", "", abstract).strip()
                if cleaned_abstract:
                    abstract = cleaned_abstract
                else:
                    print("Abstract became empty after cleaning tags.")
                    abstract = (
                        "Abstract found but empty after cleaning."  # Or keep original?
                    )
            # Check if abstract is just the default message from CrossRef function
            if abstract == "Abstract not found.":
                abstract = None  # Reset if it was the default message

        # 4. Set default message if still no abstract
        if abstract is None:
            print(f"Abstract could not be retrieved for DOI {doi} from any source.")
            article["summary"] = (
                "Abstract not found."  # Use the default message for Gemini
            )
        else:
            print(f"Successfully obtained abstract from {source}.")
            article["summary"] = abstract  # Use the retrieved abstract

        extracted_info.append(article)
        processed_count += 1
        time.sleep(1)  # Add delay between processing articles to be polite to APIs

    print(f"\nFinished processing {processed_count} articles.")

    if extracted_info:
        print("Sending messages to Discord...")
        generate_and_send_messages(extracted_info, webhook_urls)
        # Save hash only if processing was successful and new entries were found
        save_latest_entry_hash(
            get_latest_entry_hash(feedparser.parse(rss_url)), hash_file_path
        )


if __name__ == "__main__":
    main()
