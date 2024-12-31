# -*- coding: utf-8 -*-
"""Pubmed_fetch_all_data

Modified version without OpenAI translation features.
"""

import os
from Bio import Entrez, Medline
from datetime import datetime, timedelta
import gspread
from gspread_dataframe import set_with_dataframe
from google.oauth2.service_account import Credentials
import pandas as pd
from urllib.error import HTTPError
from googleapiclient.discovery import build
from tqdm import tqdm

# --- Configuration ---
SEARCH_TERM = "p-tau"  # Customize your search term
DAYS_BACK = 7  # Number of days back to search
DATE_FORMAT = "%Y/%m/%d"
MAX_RETMAX = 200

# Google Sheets Configuration
SCOPE = [
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive",
]
CREDS_FILE = "account_credential.json"  # Path to your Google service account credentials JSON file
SHEET_NAME = "Name_of_your_Google_Sheet"  # Name of your Google Sheet
WORKSHEET_NAME = SEARCH_TERM  # Name of the worksheet within the Google Sheet

# Set your email for Entrez (use environment variable)
Entrez.email = os.environ.get("ENTREZ_EMAIL")

def fetch_pubmed_ids(search_term, days_back):
    """
    Fetches PubMed IDs for articles matching the search term within the last 'days_back' days.
    """
    with tqdm(total=1, desc="Fetching PubMed IDs") as pbar:
        handle = Entrez.esearch(
            db="pubmed",
            term=search_term,
            reldate=days_back,
            datetype="pdat",
            retmax=MAX_RETMAX,
            usehistory="y"
        )
        search_results = Entrez.read(handle)
        handle.close()
        pbar.update(1)
    return search_results["IdList"], search_results["WebEnv"], search_results["QueryKey"]

def fetch_pubmed_records(id_list, webenv, query_key):
    """
    Fetches PubMed records (details) for the given list of IDs.
    """
    if not id_list:
        return []

    records = []
    with tqdm(total=len(id_list), desc="Fetching PubMed Records") as pbar:
        try:
            batch_size = 10
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

def load_to_google_sheet(df, sheet_name, worksheet_name, creds_file, scope):
    """
    Loads a Pandas DataFrame to a Google Sheet.
    """
    creds = Credentials.from_service_account_file(creds_file, scopes=scope)
    client = gspread.authorize(creds)

    try:
        sheet = client.open(sheet_name)
    except gspread.exceptions.SpreadsheetNotFound:
        print(f"Spreadsheet '{sheet_name}' not found. Creating a new one...")
        sheet = client.create(sheet_name)
        file_id = sheet.id
        folder_id = "YOUR_FOLDER_ID"  # Replace with your actual folder ID

        drive_service = build('drive', 'v3', credentials=creds)
        file = drive_service.files().get(fileId=file_id, fields='parents').execute()
        previous_parents = ",".join(file.get('parents'))
        file = drive_service.files().update(fileId=file_id,
                                          addParents=folder_id,
                                          removeParents=previous_parents,
                                          fields='id, parents').execute()

    try:
        worksheet = sheet.worksheet(worksheet_name)
    except gspread.exceptions.WorksheetNotFound:
        print(f"Worksheet '{worksheet_name}' not found. Creating a new one...")
        worksheet = sheet.add_worksheet(title=worksheet_name, rows="100", cols="20")

    with tqdm(total=3, desc="Loading to Google Sheet") as pbar:
        worksheet.clear()
        pbar.update(1)
        set_with_dataframe(worksheet, df)
        pbar.update(1)

        pubmed_links = [f'=HYPERLINK("https://pubmed.ncbi.nlm.nih.gov/{pmid}", "{pmid}")' for pmid in df['pmid']]
        worksheet.update('A2:A{}'.format(len(df) + 1), [[link] for link in pubmed_links], value_input_option='USER_ENTERED')
        pbar.update(1)

    print(f"Data successfully loaded to '{sheet_name}' -> '{worksheet_name}'")

def pubmed_data_pipeline(search_term, days_back):
    """
    Main function for the PubMed data processing pipeline.
    Fetches, processes, formats, and loads PubMed data into a Google Sheet.
    """
    id_list, webenv, query_key = fetch_pubmed_ids(search_term, days_back)
    records = fetch_pubmed_records(id_list, webenv, query_key)

    print("Formatting records...")
    with tqdm(total=len(records), desc="Formatting records") as pbar:
        formatted_records = []
        for record in records:
            formatted_records.append(format_record(record))
            pbar.update(1)

    print("Generating summary...")
    summary = generate_summary(formatted_records)

    if formatted_records:
        # Create a Pandas DataFrame from the formatted records
        df = pd.DataFrame(formatted_records)

        # Load the DataFrame to Google Sheets
        print("Loading data to Google Sheets...")
        load_to_google_sheet(
            df, SHEET_NAME, WORKSHEET_NAME, CREDS_FILE, SCOPE
        )

    return formatted_records, summary

# --- Example Usage ---
if __name__ == "__main__":
    formatted_records, summary = pubmed_data_pipeline(SEARCH_TERM, DAYS_BACK)
    print("-" * 30)
    print(summary)
    print("-" * 30)
