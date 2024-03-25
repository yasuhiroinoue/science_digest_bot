import feedparser
import requests
import hashlib
import os
import time
from Bio import Entrez
from anthropic import AnthropicVertex
from dotenv import load_dotenv
from dateutil.parser import parse as parse_date  # 日付を解析するため

# Load environment variables from a .env file
load_dotenv()

# Replace hard-coded values with values from environment variables
model = os.getenv("MODEL", "default_model_value")
region = os.getenv("REGION", "us-central1")
project_id = os.getenv("PROJECT_ID", "default_project_id")
rss_url = os.getenv("RSS_URL", "default_rss_url")
webhook_url = os.getenv("WEBHOOK_URL", "default_webhook_url")
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
        return [], False  # 新しい記事がないことを示す
    
    articles = []
    for entry in feed.entries:
        doi = entry.get('dc_identifier', '').split('doi:')[-1].strip()
        title = entry.title
        updated = entry.get('updated', 'No date available')  # updated キーから日付を取得
        parsed_date = parse_date(updated) if updated != 'No date available' else 'Date unavailable'  # 日付を解析する
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
        return "DOI not available, possibly a news item or content without a DOI."
    
    try:
        time.sleep(1/2)  # NCBIの使用制限に従う
        search_result = Entrez.esearch(db="pubmed", term=doi, retmax=10)
        record = Entrez.read(search_result)
        
        if not record["IdList"]:
            return "No corresponding PubMed entry found."
        
        id_list = record["IdList"]
        details = Entrez.efetch(db="pubmed", id=id_list[0], retmode="xml", validate=False)
        details_record = Entrez.read(details)
        
        if 'PubmedArticle' not in details_record or not details_record['PubmedArticle']:
            return "Abstract not found."
        
        article = details_record['PubmedArticle'][0]
        if 'MedlineCitation' in article and 'Article' in article['MedlineCitation'] and 'Abstract' in article['MedlineCitation']['Article'] and 'AbstractText' in article['MedlineCitation']['Article']['Abstract']:
            abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            return abstract
        else:
            return "Abstract not found."
    except Exception as e:
        return f"Error retrieving abstract: {str(e)}"

def generate_and_send_messages(extracted_info, model, region, project_id, webhook_url):
    client = AnthropicVertex(region=region, project_id=project_id)
    
    for info in extracted_info:
        content_message = f"Title: {info['Title']}\nAuthors: {info['Authors']}\nDate: {info['Date']}\nSummary: {info['summary']}\nURL: {info['URL']}"

        # Send request to Anthropics API to generate content
        message = client.messages.create(
            max_tokens=4096,
            system="""
    Imagine you are a highly skilled English-Japanese translator with a specialization in academic content. 
    Your task is to reformat the information provided below for an academic audience. Keep the 'Title' and 'Authors' 
    sections in English to maintain the original context and recognition. For the 'Summary', translate the essence 
    of the work into Japanese, succinctly capturing the main points in no more than 50 characters. This limit may 
    be adjusted if the complexity of the content demands it, but strive for conciseness. Finally, append the relevant URL 
    at the end, ensuring the journal name 'Science' remains in its original language for consistency. 
    Structure your response according to the Markdown format detailed below, ensuring each element is clearly distinguished 
    for readability and accessibility:
    
    - **Title**: {Title}
    - **Authors**: {Authors}
    - **Date**: {Date}
    - **Summary**: {summary}
    - **URL**: {url}
    """,
            messages=[
                {
                    "role": "user",
                    "content": content_message,
                }
            ],
            model=model,
        )
        
        # Assume message.content[0].text contains the generated message for Discord
        discord_message = message.content[0].text
        
        # Post the generated content to Discord
        payload = {"content": discord_message}
        response = requests.post(webhook_url, json=payload)
        
        if response.status_code not in [200, 204]:
            print(f"Error: {response.status_code} - {response.reason}")

        # Respect rate limits
        time.sleep(2)  # Adjust based on your rate limit policies
        
def main():
    articles, has_new_entries = get_latest_articles_from_rss(rss_url, hash_file_path)
    if not has_new_entries:
        return
    
    extracted_info = []
    for article in articles:  
        abstract = get_abstract_from_pubmed(article['DOI'])
        article['summary'] = abstract
        extracted_info.append(article)

    if extracted_info:
        generate_and_send_messages(extracted_info, model, region, project_id, webhook_url)
        save_latest_entry_hash(get_latest_entry_hash(feedparser.parse(rss_url)), hash_file_path)

# Uncomment for execution in a real environment
if __name__ == "__main__":
    main()
