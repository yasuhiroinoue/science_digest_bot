# Science Article Digest Bot

This bot automatically fetches the latest articles from the *Science* journal RSS feed, retrieves their abstracts from Crossref, summarizes the abstracts in short using the Claude3.5-haiku model, and posts these summaries to a specified Discord channel.

## Setup

### Prerequisites

- Python 3.8 or newer
- An active internet connection
- An access to VertexAI (for Claude3)
- A Discord server with a webhook URL

### Installation

1. **Clone the Repository:**

```bash
git clone <repository-url>
cd <repository-name>
```

2. **Install Dependencies:**

Create a virtual environment and install the required Python packages.

```bash
python -m venv venv
source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
pip install -r requirements.txt
```

3. **Configure Environment Variables:**

Copy the `.env.example` file to a new file named `.env`, and fill in the values:

```plaintext
MODEL=<your-anthropic-model-identifier>
REGION=<your-anthropic-region>
PROJECT_ID=<your-anthropic-project-id>
RSS_URL=<science-journal-rss-feed-url>
WEBHOOK_URL=<your-discord-webhook-url>
HASH_FILE_PATH=./latest_entry_hash.txt
```

Ensure you replace `<placeholder>` with your actual data.

## Usage

Execute the script to fetch the latest articles, summarize them, and post to Discord:

```bash
python science.py
```

### Cron Job

For automatic daily updates, consider adding a cron job:

```bash
0 9 * * 5 /path/to/your/venv/python /path/to/science.py
```

This example runs the script every Friday at 9:00 AM.

## How It Works

1. **Fetch Articles:** Uses RSS to get the latest articles from *Science*.
2. **Retrieve Abstracts:** Looks up articles by DOI in PubMed to get abstracts.
3. **Summarize with AI:** Sends abstracts to Claude3 to generate short summaries.
4. **Post to Discord:** Formats the summaries and posts them to a Discord channel.

### Note

The script respects API rate limits by including pauses between requests. Adjust these as necessary based on your API usage policies.

## Contributing

Contributions to improve the script or documentation are welcome. Please follow the standard GitHub pull request process to propose changes.
