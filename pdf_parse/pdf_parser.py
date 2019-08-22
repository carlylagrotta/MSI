#!usr/bin/python3
import PyPDF2 as pdf
import tabula as tb

############# PDF Parser Class #############
'''
Written by Dan Lee, 2019 Burke Lab
- Extract data from PDF
   - Identify & Read tables to specified format
   - Identify graphs and approximate data from visuals
'''

class PDF_Parser(object):
    def __init__(self, path:str = ''):
        self.fpath      = path
        f               = open(path, 'rb')
        self.doc_reader = pdf.PdfFileReader(f)
        self.meta       = self.doc_reader.getDocumentInfo()
        self.num_pages  = self.doc_reader.getNumPages()

    # read_table (num= page_number, path= filename to store dataframe)
    def read_table(self, num:int, path:str = ''):
        page = self.doc_reader.getPage(page)
        if Table in 


