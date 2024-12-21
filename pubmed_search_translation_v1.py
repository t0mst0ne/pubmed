# -*- coding: utf-8 -*-

# If in mac pip install gspread gspread-dataframe google-api-python-client google-auth-oauthlib google-auth-httplib2

import os
from Bio import Entrez, Medline
from datetime import datetime, timedelta
import gspread
from gspread_dataframe import set_with_dataframe
from google.oauth2.service_account import Credentials
import pandas as pd
from urllib.error import HTTPError
from googleapiclient.discovery import build
from openai import OpenAI
from tqdm import tqdm

# Define the topic for which to generate PubMed search phrases
search_topic = "dementia treatment"  # **CHANGE THIS** to your desired topic
days_back = 28
max_retmax = 200



# Google Sheets Configuration
SCOPE = [
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive",
]
CREDS_FILE = "sublime-state-xxx.json"  # **CHANGE THIS** Path to your Google service account credentials JSON file
SHEET_NAME = "Pubmed_search_Dementia"  # **CHANGE THIS** Name of your Google Sheet
WORKSHEET_NAME = "dementia_treatment"  # **CHANGE THIS** Name of the worksheet within the Google Sheet

# Set your email for Entrez (use environment variable)
Entrez.email = os.environ.get("ENTREZ_EMAIL")  # **MAKE SURE TO SET THIS ENVIRONMENT VARIABLE**

# OpenAI API Key
OPENAI_API_KEY = 'openAIkey'


# --- Functions for Core Processing ---

def generate_pubmed_search_phrases(topic, num_phrases=3):
    """
    Generates PubMed search phrases related to a given topic using the OpenAI API.

    Args:
        topic (str): The topic for which to generate search phrases.
        num_phrases (int): The number of search phrases to generate.

    Returns:
        list: A list of generated PubMed search phrases.
    """
    client = OpenAI(api_key=OPENAI_API_KEY)
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "You are a helpful assistant that generates PubMed search phrases."},
                {"role": "user", "content": f"Generate {num_phrases} PubMed search phrases related to the topic: '{topic}'. Each search phrase should be optimized for finding relevant research articles on PubMed. Please provide only the search phrases in your response, with each phrase on a new line."}
            ]
        )
        phrases_raw = response.choices[0].message.content.strip()
        phrases = phrases_raw.split('\n')  # Split into individual phrases
        print (phrases)
        return phrases
    except Exception as e:
        print(f"Error generating PubMed search phrases: {e}")
        return []

def fetch_pubmed_ids(search_term, days_back, max_retmax=10):
    """
    Fetches PubMed IDs for articles matching the search term within the last 'days_back' days.

    Args:
        search_term (str): The search term for PubMed.
        days_back (int): Number of days back to search.
        max_retmax (int): Maximum number of results to retrieve.

    Returns:
        tuple: A tuple containing a list of PubMed IDs, WebEnv, and QueryKey.
    """
    with tqdm(total=1, desc=f"Fetching PubMed IDs for '{search_term}'") as pbar:
        try:
            handle = Entrez.esearch(
                db="pubmed",
                term=search_term,
                reldate=days_back,
                datetype="pdat",
                retmax=max_retmax,
                usehistory="y"
            )
            search_results = Entrez.read(handle)
            handle.close()
            pbar.update(1)
            return search_results["IdList"], search_results["WebEnv"], search_results["QueryKey"]
        except HTTPError as e:
            print(f"Error fetching PubMed IDs for '{search_term}': {e}")
            return [], None, None

def fetch_pubmed_records(id_list, webenv, query_key):
    """
    Fetches PubMed records (details) for the given list of IDs.

    Args:
        id_list (list): A list of PubMed IDs.
        webenv (str): Web environment key for Entrez history server.
        query_key (str): Query key for Entrez history server.

    Returns:
        list: A list of Medline records (dictionaries).
    """
    if not id_list:
        return []

    records = []
    with tqdm(total=len(id_list), desc="Fetching PubMed Records") as pbar:
        try:
            batch_size = 10  # Fetch in batches
            for i in range(0, len(id_list), batch_size):
                batch_ids = id_list[i:i + batch_size]
                handle = Entrez.efetch(
                    db="pubmed",
                    id=batch_ids,
                    rettype="medline",
                    retmode="text",
                    webenv=webenv,
                    query_key=query_key,
                )
                batch_records = Medline.parse(handle)
                for rec in batch_records:
                    records.append(rec)
                    pbar.update(1)
                handle.close()
        except HTTPError as e:
            print(f"Error fetching records: {e}")
            if e.code == 400:
                print("The WebEnv or QueryKey might be invalid or expired.")
            return []

    return records

def format_record(record):
    """
    Formats a single Medline record into a more readable format.

    Args:
        record (dict): A Medline record.

    Returns:
        dict: A formatted record.
    """
    formatted_record = {
        "pmid": record.get("PMID", ""),
        "title": record.get("TI", ""),
        "abstract": record.get("AB", ""),
        "authors": ", ".join(record.get("AU", [])),
        "journal": record.get("JT", ""),
        "pubdate": record.get("DP", ""),
    }
    return formatted_record

def generate_summary(formatted_records):
    """
    Generates a summary from a list of formatted records.

    Args:
        formatted_records (list): A list of formatted records.

    Returns:
        str: A summary string.
    """
    summary = ""
    if not formatted_records:
        summary = "No new articles found for the given search criteria today."
        return summary

    for record in formatted_records:
        summary += f"**{record['title']}**\n"
        summary += f"PMID: {record['pmid']}\n"
        summary += f"Authors: {record['authors']}\n"
        summary += f"Journal: {record['journal']}\n"
        summary += f"Publication Date: {record['pubdate']}\n"
        if record["abstract"]:
            summary += f"Abstract: {record['abstract']}\n\n"
        else:
            summary += "Abstract: N/A\n\n"
    return summary

def translate_to_traditional_chinese(text):
    """
    Translates text to Traditional Chinese using the OpenAI API.

    Args:
        text (str): The text to translate.

    Returns:
        str: The translated text, or None if an error occurred.
    """
    client = OpenAI(api_key=OPENAI_API_KEY)
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "You are a helpful assistant that translates English to Traditional Chinese as used in Taiwan."},
                {"role": "user", "content": f"Translate the following text to Traditional Chinese (Taiwan):\n\n{text}"}
            ]
        )
        translated_text = response.choices[0].message.content.strip()
        return translated_text
    except Exception as e:
        print(f"Error during translation: {e}")
        return None

def load_to_google_sheet(df, sheet_name, worksheet_name, creds_file, scope):
    """
    Loads a Pandas DataFrame to a Google Sheet.

    Args:
        df (pd.DataFrame): The DataFrame to load.
        sheet_name (str): The name of the Google Sheet.
        worksheet_name (str): The name of the worksheet.
        creds_file (str): Path to the credentials JSON file.
        scope (list): List of authorization scopes.
    """
    creds = Credentials.from_service_account_file(creds_file, scopes=scope)
    client = gspread.authorize(creds)

    # Open the Google Sheet
    try:
        sheet = client.open(sheet_name)
    except gspread.exceptions.SpreadsheetNotFound:
        print(f"Spreadsheet '{sheet_name}' not found. Creating a new one...")
        sheet = client.create(sheet_name)
        # Get the file ID of the newly created spreadsheet
        file_id = sheet.id
        # Get the folder ID where you want to move the file
        folder_id = "YOUR_FOLDER_ID"  # Replace with your actual folder ID

        # Move the file using the Drive API
        drive_service = build('drive', 'v3', credentials=creds)
        file = drive_service.files().get(fileId=file_id, fields='parents').execute()
        previous_parents = ",".join(file.get('parents'))
        file = drive_service.files().update(fileId=file_id,
                                            addParents=folder_id,
                                            removeParents=previous_parents,
                                            fields='id, parents').execute()


    # Select the worksheet
    try:
        worksheet = sheet.worksheet(worksheet_name)
    except gspread.exceptions.WorksheetNotFound:
        print(f"Worksheet '{worksheet_name}' not found. Creating a new one...")
        worksheet = sheet.add_worksheet(title=worksheet_name, rows="100", cols="20")  # Adjust rows/cols as needed

    with tqdm(total=3, desc="Loading to Google Sheet") as pbar:
        # Clear existing data and update with DataFrame
        worksheet.clear()
        pbar.update(1)
        set_with_dataframe(worksheet, df)
        pbar.update(1)

        # Add hyperlinks to PMID column
        pubmed_links = [f'=HYPERLINK("https://pubmed.ncbi.nlm.nih.gov/{pmid}", "{pmid}")' for pmid in df['pmid']]
        worksheet.update('A2:A{}'.format(len(df) + 1), [[link] for link in pubmed_links], value_input_option='USER_ENTERED')
        pbar.update(1)

    print(f"Data successfully loaded to '{sheet_name}' -> '{worksheet_name}'")

def pubmed_data_pipeline(topic, days_back, max_retmax=10):
    """
    Main function for the PubMed data processing pipeline.
    Generates search phrases, fetches data for each phrase, combines results,
    and loads data into a Google Sheet.
    """
    generated_phrases = generate_pubmed_search_phrases(topic)

    all_formatted_records = []

    for phrase in generated_phrases:
        id_list, webenv, query_key = fetch_pubmed_ids(phrase, days_back, max_retmax)
        records = fetch_pubmed_records(id_list, webenv, query_key)

        print(f"Formatting records for phrase: {phrase}...")
        with tqdm(total=len(records), desc=f"Formatting records for '{phrase}'") as pbar:
            formatted_records = []
            for record in records:
                formatted_records.append(format_record(record))
                pbar.update(1)

            all_formatted_records.extend(formatted_records)

    print("Generating summary...")
    summary = generate_summary(all_formatted_records)

    if all_formatted_records:
        df = pd.DataFrame(all_formatted_records)

        # Remove duplicates based on PMID
        df.drop_duplicates(subset='pmid', keep='first', inplace=True)

        print("Translating abstracts to Traditional Chinese...")
        with tqdm(total=len(df), desc="Translating abstracts") as pbar:
            df["abstract_zh_TW"] = df["abstract"].apply(lambda x: translate_to_traditional_chinese(x) if x else None)
            pbar.update(len(df))

        print("Loading data to Google Sheets...")
        load_to_google_sheet(df, SHEET_NAME, WORKSHEET_NAME, CREDS_FILE, SCOPE)

    return all_formatted_records, summary

# --- Example Usage ---
if __name__ == "__main__":
    # Set OpenAI API key
    OpenAI.api_key = OPENAI_API_KEY

    formatted_records, summary = pubmed_data_pipeline(search_topic, days_back, max_retmax)
    print("-" * 30)
    print(summary)
    print("-" * 30)

