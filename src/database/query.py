import sqlite3

conn = sqlite3.connect("MARK.db")
cursor = conn.cursor()


# ADD CODE HERE TO parse FITS 
# image headers and automatically store the data in this database
# follow format of our schematic




conn.close()