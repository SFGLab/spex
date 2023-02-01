import gspread
import datetime
from oauth2client.service_account import ServiceAccountCredentials
import traceback
import time
import math
import random
import logging

class GSheetLogger:
    def __init__(self, google_sheets_credentials_filename, spreadsheet_id):
        scope = ['https://spreadsheets.google.com/feeds',
                 'https://www.googleapis.com/auth/drive']
        self.credentials = ServiceAccountCredentials.from_json_keyfile_name(google_sheets_credentials_filename, scope)
        self.spreadsheet_id = spreadsheet_id

    def log(self, worksheet, params, depth=0):
        try:
            if(depth > 0):
                time.sleep(math.pow(2, depth)+random.random())
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
            print("Logged")

        except Exception:
            if(depth > 7):                
                logging.info("Exception while adding the row to Google Sheet")
                print(traceback.format_exc())
            else:
                logging.info("Exception while adding the row to Google Sheet - retry #" + str(depth+1) + "/8")
                self.log(worksheet, params, depth+1)
