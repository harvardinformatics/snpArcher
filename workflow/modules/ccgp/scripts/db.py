import pymongo
from dotenv import load_dotenv
import os


def get_mongo_client() -> pymongo.MongoClient:
    """Utility function to get MongoDB client."""
    load_dotenv()
    DB_CONNECT_STRING = os.getenv("DB_CONNECT_STRING")
    client = pymongo.MongoClient(DB_CONNECT_STRING)
    return client
