import requests
import xml.etree.ElementTree as ET
from time import sleep
import pandas as pd

def search_pubmed(term, retmax=100):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    search_url = f"{base_url}esearch.fcgi?db=pubmed&term={term}&retmax={retmax}&usehistory=y"
    
    response = requests.get(search_url)
    root = ET.fromstring(response.content)
    
    webenv = root.find("WebEnv").text
    query_key = root.find("QueryKey").text
    count = int(root.find("Count").text)
    
    return webenv, query_key, count

def fetch_articles(webenv, query_key, retmax=100, retstart=0):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    fetch_url = f"{base_url}efetch.fcgi?db=pubmed&WebEnv={webenv}&query_key={query_key}&retmax={retmax}&retstart={retstart}&retmode=xml"
    
    response = requests.get(fetch_url)
    root = ET.fromstring(response.content)
    
    articles = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.find(".//PMID").text
        title = article.find(".//ArticleTitle").text
        journal = article.find(".//Journal/Title").text
        
        # Extract publication year
        pub_date = article.find(".//PubDate")
        year = pub_date.find("Year")
        pub_year = year.text if year is not None else "N/A"
        
        # Extract DOI
        article_id_list = article.find(".//ArticleIdList")
        doi = "N/A"
        if article_id_list is not None:
            for article_id in article_id_list.findall("ArticleId"):
                if article_id.get("IdType") == "doi":
                    doi = article_id.text
                    break
        
        articles.append({
            "pmid": pmid,
            "title": title,
            "journal": journal,
            "year": pub_year,
            "doi": doi
        })
    
    return articles

def main():
    search_term = '"fall+detection"'
    webenv, query_key, total_count = search_pubmed(search_term)
    
    all_articles = []
    retmax = 100
    
    for retstart in range(0, total_count, retmax):
        articles = fetch_articles(webenv, query_key, retmax, retstart)
        all_articles.extend(articles)
        print(f"Retrieved {len(all_articles)} articles out of {total_count}")
        sleep(1)  # To comply with NCBI's usage guidelines
    
    # Create pandas DataFrame
    df = pd.DataFrame(all_articles)
    
    # Display the first few rows of the DataFrame
    print(df.head())
    
    # Optionally, save the DataFrame to a CSV file
    df.to_csv("fall_detection_articles.csv", index=False)
    print("DataFrame saved to 'fall_detection_articles.csv'")

if __name__ == "__main__":
    main()
