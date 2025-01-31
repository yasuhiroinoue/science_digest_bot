import feedparser
import requests
import hashlib
import os
import time
from dotenv import load_dotenv
from dateutil.parser import parse as parse_date  # To parse date
from google import genai
from google.genai import types


# Load environment variables from a .env file
load_dotenv()

# Replace hard-coded values with values from environment variables
region = os.getenv("REGION", "us-central1")
project_id = os.getenv("PROJECT_ID", "default_project_id")
rss_url = os.getenv("RSS_URL", "default_rss_url")
webhook_urls = os.environ.get("WEBHOOK_URLS").split(",")
hash_file_path = os.getenv("HASH_FILE_PATH", "latest_entry_hash.txt")


# Gemini
MODEL_ID = "gemini-2.0-flash-exp"

generate_content_config = types.GenerateContentConfig(
    temperature = 1,
    top_p = 0.95,
    # top_k= 64,
    max_output_tokens = 8192,
    system_instruction="""
あなたは高度な英語-日本語翻訳者です。提供された英語の論文要旨を学術的な聴衆に向けて日本語に要約してください。要約は100文字以内に抑え、論文の核心を簡潔に表現してください。文字数の制限は論文の内容の複雑さによって柔軟に調整可能です。要約が不可能な場合、または「Abstract not found.」と入力された場合は、「要約なし」と返してください。論文要旨が短い場合、直接的な日本語訳を提供してください。プロセスの結果として、要約または翻訳されたテキストのみを返すようにしてください。「要旨:」という見出しを付けることは禁止します。
""",
)

gemini_client = genai.Client(
    vertexai=True, project=project_id, location=region
)

# Gemini response
def generate_text(prompt):
    response = gemini_client.models.generate_content(
    model=MODEL_ID, contents=prompt, config=generate_content_config
    )
    return response.text

def get_latest_entry_hash(rss_data):
    latest_entry = rss_data['entries'][0]
    unique_str = latest_entry.get('title', '') + latest_entry.get('published', '')
    return hashlib.md5(unique_str.encode()).hexdigest()


def save_latest_entry_hash(hash_str, hash_file_path):
    with open(hash_file_path, 'w') as file:
        file.write(hash_str)


def read_latest_entry_hash(hash_file_path):
    if os.path.exists(hash_file_path):
        with open(hash_file_path, 'r') as file:
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
        doi = entry.get('dc_identifier', '').split('doi:')[-1].strip()
        title = entry.title
        updated = entry.get('updated', 'No date available')  # Get date from the updated key
        parsed_date = parse_date(updated) if updated != 'No date available' else 'Date unavailable'  # parse date
        authors = entry.get('author', 'Authors not listed')
        url = entry.link
        articles.append({
            'DOI': doi,
            'Title': title,
            'Authors': authors,
            'Date': parsed_date,  # 記事の辞書に日付を追加
            'URL': url
        })
    return articles, True


def get_abstract_from_crossref(doi):
    if not doi:
        return "Abstract not found."
    
    try:
        time.sleep(1)  # API制限を考慮してスリープを入れる
        api_url = f"https://api.crossref.org/works/{doi}"
        headers = {
            'User-Agent': 'DOIScraperBot/1.0 (mailto:your-email@example.com)',
            'Accept': 'application/json'
        }
        response = requests.get(api_url, headers=headers)
        response.raise_for_status()
        data = response.json()

        if 'message' not in data:
            return "Unexpected API response format"

        message = data['message']
        abstract = message.get('abstract', "Abstract not found.")
        # For debug
        # print(abstract)
        return abstract
    except requests.RequestException as e:
        return f"Error retrieving abstract: {str(e)}"


def format_authors(info):
    authors = info['Authors'].replace('\n', '')  # Remove any newline characters
    author_list = []
    
    # Split the author names into a list
    author_names = authors.split(", ")
    
    if len(author_names) > 4:
        # Select first three authors and append 'et al.'
        formatted_authors = f"{author_names[0]}, {author_names[1]}, {author_names[2]} et al."
    else:
        # Join all authors because there are three or less
        formatted_authors = ", ".join(author_names)
    
    # Append the formatted string to author_list
    author_list.append(formatted_authors)
    
    return author_list


def format_to_markdown(content_message):
    # Extract the information from the content message
    info = {}
    for line in content_message.split('\n'):
        if line.strip():  # Check if line is not empty
            key, value = line.split(':', 1)
            info[key.strip()] = value.strip()
    
    # Remove blank lines from the Summary
    if 'Summary' in info:
        lines = info['Summary'].split('\n')
        filtered_summary = '\n'.join([line for line in lines if line.strip()])
        info['Summary'] = filtered_summary

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
        retry_delay = 5  # 再試行の間隔（秒）

        for url in webhook_urls:
            for attempt in range(max_retries):
                response = requests.post(url, json=payload)
                
                if response.status_code in [200, 204]:
                    break  # 成功したらループを抜ける
                else:
                    print(f"Error: {response.status_code} - {response.reason}. Retrying...")
                    print(discord_message)
                    time.sleep(retry_delay)  # 指定された時間待つ
            else:
                # すべての再試行が失敗した場合
                print(f"Failed to post after {max_retries} attempts.")
        
        time.sleep(30)  # Adjust based on your rate limit policies


def main():
    articles, has_new_entries = get_latest_articles_from_rss(rss_url, hash_file_path)
    if not has_new_entries:
        return
    
    extracted_info = []
    
    for article in articles[:]:
        abstract = get_abstract_from_crossref(article['DOI'])
        article['summary'] = abstract
        extracted_info.append(article)

    if extracted_info:
        generate_and_send_messages(extracted_info, webhook_urls)
        save_latest_entry_hash(get_latest_entry_hash(feedparser.parse(rss_url)), hash_file_path)


if __name__ == "__main__":
    main()
