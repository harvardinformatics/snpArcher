from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
import google.auth
import io
from googleapiclient.http import MediaIoBaseDownload, MediaFileUpload
from pathlib import Path
from dotenv import load_dotenv

SCOPES = ["https://www.googleapis.com/auth/drive"]


class CCGPDrive:
    """Class for interacting with CCGP Data Wrangling drive."""

    def __init__(self) -> None:
        """Gets credentials and builds google drive service."""
        load_dotenv()
        self.creds, _ = google.auth.default(scopes=SCOPES)
        self.service = build("drive", "v3", credentials=self.creds)

    def _get_files_list_response(self, query: str) -> list[dict]:
        """ "Helper function that sends creates and send query to Drive API and returns response."""
        page_token = None
        result = []
        while True:
            response = (
                self.service.files()
                .list(
                    q=query,
                    spaces="drive",
                    fields="nextPageToken, files(id, name, modifiedTime, mimeType)",
                    pageToken=page_token,
                )
                .execute()
            )
            files = response.get("files", [])
            for file in files:
                result.append(file)

            page_token = response.get(
                "nextPageToken", None
            )  # Drive API iterates each page of the drive, so we have to do this to make sure we search each page. Doesn't matter for this case but its best practices
            if page_token is None:
                break
        return result

    def get_folder_id(self, folder: str) -> str:
        query = f"name = '{folder}' and mimeType = 'application/vnd.google-apps.folder'"
        found = self._get_files_list_response(query)
        if len(found) == 0:
            return False
        result = found[0].get("id")
        return result

    def list_files_from_folder(self, folder: str) -> list[dict]:
        folder_id = self.get_folder_id(folder)
        query = f"'{folder_id}' in parents"
        found = self._get_files_list_response(query)
        return found

    def download_files(self, *files: dict) -> None:
        """Downloads files to current directory."""
        for file in files:
            if Path(file["name"]).exists():
                print(
                    "File: "
                    + "'"
                    + file["name"]
                    + "'"
                    + " already exists, skipping download."
                )
                continue
            if file.get("mimeType") == "application/vnd.google-apps.spreadsheet":
                request = self.service.files().export_media(
                    fileId=file.get("id"), mimeType="text/tab-separated-values"
                )
                file["name"]
            else:
                request = self.service.files().get_media(fileId=file.get("id"))
            fh = io.FileIO(file.get("name"), "wb")
            downloader = MediaIoBaseDownload(fh, request)
            done = False
            print("Downloading file " + "'" + file["name"] + "," + file["id"] + "'")
            while done is False:
                status, done = downloader.next_chunk()

    def upload_file(
        self, file: Path, folder_name: str = None, folder_id: str = None
    ) -> None:
        if folder_name is not None and folder_id is None:
            folder_id = self.get_folder_id(u)
            file_metadata = {"name": file.name, "parents": [folder_id]}
            if folder_id == False:
                folder_id = self.create_folder(folder_name)
                file_metadata = {"name": file.name, "parents": [folder_id]}
        elif folder_id is not None:
            file_metadata = {"name": file.name, "parents": [folder_id]}
        else:
            file_metadata = {"name": file.name}
        media = MediaFileUpload(file)
        up = (
            self.service.files()
            .create(body=file_metadata, media_body=media, fields="id")
            .execute()
        )
        print("Uploaded file: " + "'" + file.name + "'")

    def create_folder(self, folder_name: str, parent_id: str = None) -> str:
        """ "Creates folder and returns id"""

        if parent_id is None:
            parent_id = self.get_folder_id("Project Results")

        file_metadata = {
            "name": folder_name,
            "parents": [parent_id],
            "mimeType": "application/vnd.google-apps.folder",
        }
        file = (
            self.service.files()
            .create(body=file_metadata, fields="id", supportsAllDrives=True)
            .execute()
        )
        return file.get("id")


def main():
    g = CCGPDrive()
    hi = g.create_folder("hi")
    bye = g.create_folder("bye", hi)


if __name__ == "__main__":
    main()
