import feedparser
import requests
import hashlib
import os
import time
from Bio import Entrez
from anthropic import AnthropicVertex
from dotenv import load_dotenv
from dateutil.parser import parse as parse_date  # To parse date

# Load environment variables from a .env file
load_dotenv()

# Replace hard-coded values with values from environment variables
LLM_model = os.getenv("MODEL", "claude-3-haiku@20240307")
region = os.getenv("REGION", "us-central1")
project_id = os.getenv("PROJECT_ID", "default_project_id")
rss_url = os.getenv("RSS_URL", "default_rss_url")
# webhook_url = os.getenv("WEBHOOK_URL", "default_webhook_url")
webhook_urls = os.environ.get("WEBHOOK_URLS").split(",")
hash_file_path = os.getenv("HASH_FILE_PATH", "latest_entry_hash.txt")
Entrez.email = os.getenv("ENTREZ_EMAIL", "default_email@gmail.com")

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

def get_abstract_from_pubmed(doi):
    if not doi:
        return "Abstract not found."
    
    try:
        time.sleep(1)  # NCBIの使用制限に従う
        # search_result = Entrez.esearch(db="pubmed", term=doi, retmax=1)
        search_result = Entrez.esearch(db="pubmed", term=f"{doi}[doi]", retmax=10)
        record = Entrez.read(search_result)
        
        if not record["IdList"]:
            return "No corresponding PubMed entry found."
        
        id_list = record["IdList"]
        details = Entrez.efetch(db="pubmed", id=id_list[0], retmode="xml", validate=False)
        details_record = Entrez.read(details)
        
        if 'PubmedArticle' not in details_record or not details_record['PubmedArticle']:
            return "Abstract not found."
        
        article = details_record['PubmedArticle'][0]
        # For debug
        # print(article)
        if 'MedlineCitation' in article and 'Article' in article['MedlineCitation'] and 'Abstract' in article['MedlineCitation']['Article'] and 'AbstractText' in article['MedlineCitation']['Article']['Abstract']:
            abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            print(abstract)
            return abstract
        else:
            return "Abstract not found."
    except Exception as e:
        return f"Error retrieving abstract: {str(e)}"
    
def format_authors(info):
    authors = info['Authors'].replace('\n', '')  # Remove any newline characters
    # authors = info['Authors']
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

# def format_to_markdown(content_message):
#     # Extract the information from the content message
#     info = {}
#     for line in content_message.split('\n'):
#         if line.strip():  # Check if line is not empty
#             key, value = line.split(':', 1)
#             info[key.strip()] = value.strip()
    
#     # Format the extracted information into markdown
#     markdown_format = (
#         f"- **{info['Title']}**\n"
#         f" {info['Authors']}\n"
#         f" {info['Date']}\n"
#         f" ```{info['Summary']}```\n"
#         f" {info['URL']}\n"
#         # f"- **Title**: {info['Title']}\n"
#         # f"- **Authors**: {info['Authors']}\n"
#         # f"- **Date**: {info['Date']}\n"
#         # f"- **Summary**: {info['Summary']}\n"
#         # f"- **URL**: {info['URL']}"
#     )
    
#     return markdown_format
def generate_and_send_messages(extracted_info, LLM_model, region, project_id, webhook_urls):
    client = AnthropicVertex(region=region, project_id=project_id)
    
    count = 1
    for info in extracted_info:
        print(count)
        count = count + 1
        # content_message = f"Title: {info['Title']}\nAuthors: {info['Authors']}\nDate: {info['Date']}\nSummary: {info['summary']}\nURL: {info['URL']}"
        content_message = f"Summary: {info['summary']}\n"

        # Send request to Anthropics API to generate content
        message = client.messages.create(
            max_tokens=2048,
            # top_p = 0.95,
            # temperature=0,
            system="""
あなたは高度な英語-日本語翻訳者です。提供された英語の論文要旨を学術的な聴衆に向けて日本語に要約してください。要約は60文字以内に抑え、論文の核心を簡潔に表現してください。文字数の制限は論文の内容の複雑さによって柔軟に調整可能です。要約が不可能な場合、または「Abstract not found.」と入力された場合は、「要約なし」と返してください。論文要旨が短い場合、直接的な日本語訳を提供してください。プロセスの結果として、要約または翻訳されたテキストのみを返すようにしてください。「要旨:」という見出しを付けることは禁止します。
""",
#             system="""
# 高度な英語-日本語翻訳者として、入力された論文要旨を学術的な聴衆向けに要約してください。要約では、論文要旨の本質を70文字以内の日本語で簡潔に伝えてください。この制限は、内容の複雑さに応じて調整可能ですが、簡潔さを心がけてください。例外処理として、"Abstract not found."のときは、「要約なし」を表示してください。論文要旨が短い場合は、その日本語訳を表示してください。プロセスの結果として要約または翻訳されたテキストのみを返してください。「要旨:」という見出しを付けてはいけません。
# """,
#             system="""
# Given an academic abstract input in English, summarize it concisely in Japanese for an academic audience. The summary should capture the essence of the abstract in no more than 70 characters, adjusted as necessary for complexity but striving for brevity. In cases where the abstract is not found, display "要約なし". If the abstract is already short, provide a direct translation in Japanese. Return only the summarized or translated text.
# """,
            messages=[
                {
                    "role": "user",
                    "content": content_message,
                }
            ],
            model=LLM_model,
        )

        abstract = message.content[0].text
        authors = format_authors(info)
        content = f"Title: {info['Title']}\nAuthors: {authors[0]}\nDate: {info['Date']}\nSummary: {abstract}\nURL: {info['URL']}"        

        print(content)
        # For debug
        # print({info['summary']})
        
        # Assume message.content[0].text contains the generated message for Discord
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

        # for url in webhook_urls:
        #     response = requests.post(url, json=payload)
        
        #     if response.status_code not in [200, 204]:
        #         print(f"Error: {response.status_code} - {response.reason}")
        
        # Respect rate limits
        time.sleep(5)  # Adjust based on your rate limit policies
        
def main():
    articles, has_new_entries = get_latest_articles_from_rss(rss_url, hash_file_path)
    if not has_new_entries:
        return
    
    extracted_info = []
    
    for article in articles[:]:
        abstract = get_abstract_from_pubmed(article['DOI'])
        article['summary'] = abstract
        extracted_info.append(article)
        # print(article['DOI'])
        # print(abstract)

    if extracted_info:
        generate_and_send_messages(extracted_info, LLM_model, region, project_id, webhook_urls)
        save_latest_entry_hash(get_latest_entry_hash(feedparser.parse(rss_url)), hash_file_path)

# Uncomment for execution in a real environment
if __name__ == "__main__":
    main()
