import gspread
import datetime
from oauth2client.service_account import ServiceAccountCredentials
import traceback


class GSheetLogger:
    def __init__(self, google_sheets_credentials_filename, spreadsheet_id):
        scope = ['https://spreadsheets.google.com/feeds',
                 'https://www.googleapis.com/auth/drive']
        self.credentials = ServiceAccountCredentials.from_json_keyfile_name(google_sheets_credentials_filename, scope)
        self.spreadsheet_id = spreadsheet_id

    def log(self, worksheet, params):
        try:
            gc = gspread.authorize(self.credentials)
            sheet = gc.open_by_key(self.spreadsheet_id).worksheet(worksheet)
            sheet_keys = sheet.row_values(1)

            params = {k.lower(): v for k, v in params.items()}
            values = list()
            for key in sheet_keys:
                key = key.lower()
                if key == "date":
                    val = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
                elif key in params:
                    val = params[key]
                else:
                    val = "?"
                if type(val) == list:
                    values.append(", ".join(val))
                else:
                    values.append(val)
            sheet.append_row(values, insert_data_option="INSERT_ROWS")

        except Exception:
            print("Exception while adding the row to Google Sheet")
            print(traceback.format_exc())
