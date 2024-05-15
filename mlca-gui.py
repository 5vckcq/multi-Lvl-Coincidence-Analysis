# to do:
# 1) add possiblity to enter constraints
# 2) show coincidence table when importing data from csv
# 3) show asf when importing data from QCA/CNA
# 4) improve rudimentary zoom function

#!/usr/bin/env python3

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import PhotoImage

import os
import math
import fitz # pdf operations from pymupdf package  -> new dependency pymupdf
from PIL import Image, ImageTk

import string # provides list of letters
import csv
import time
from datetime import datetime

import mLCA



version = "1.0"
default_num_conf = 20
default_num_var = 5
entry_list = [] # list to keep track of all generated entry widgets for table cells

pdf_width = 560
pdf_height = 640


class PDFObject:
    def __init__(self, filepath):
        # creating the file path
        self.filepath = filepath
        # opening the pdf document
        self.pdf = fitz.open(self.filepath)
        # loading the first page of the pdf document
        self.first_page = self.pdf.load_page(0)
        # getting the height and width of the first page
        self.width, self.height = self.first_page.rect.width, self.first_page.rect.height
        # zooming the page
        self.zoom = min(int(math.floor(pdf_width / self.width)), int(math.floor(pdf_height / self.height)))

    # the function for getting the page
    def get_page(self, page_num, zoom_modifier=0):
        # loading the page
        page = self.pdf.load_page(page_num)
        mat = fitz.Matrix(int(self.zoom + zoom_modifier), int(self.zoom + zoom_modifier))
        # gets the image of the page
        pix = page.get_pixmap(matrix=mat)
        
        # a variable that holds a transparent image
        px1 = pix
        # converting the image to bytes
        imgdata = px1.tobytes("ppm")
        # returning the image data
        return PhotoImage(data=imgdata)

class NotebookGridApp:
    def __init__(self, root):
        self.root = root
        self.root.title("mLCA " + version)
        
        self.current_page = 0
        self.numPages = None
        self.fileisopen = None
        self.pdf_width = pdf_width + 20
        self.pdf_height = pdf_height + 20

        # create notebook and add csv and cna tab
        self.create_notebook()       
        
        self.create_constraints()
        
        self.create_run()
        
        # add a table to the right of the notebooks
        self.num_rows = 4
        self.num_columns = 3
        self.create_outer_table()

        self.create_pdf_display()

    def update(self,new=False):
        rows = int(self.num_configurations.get()) + 1 # one additional row for level/causal ordering information
        cols = int(self.num_variables.get())
        
        if cols > 26:
            cols = 26
            print('The maximum number of variables is 26.')
            
        factor_list = list(string.ascii_uppercase)[:cols] # define factor_list as the first letters of the alphabet


        if not new:
            # delete all data
            self.table.delete(*self.table.get_children())
            
            # destroy generated entry widgets
            for widget in entry_list:
                widget.destroy()    

        # new column IDs and displayed columns
        cols = cols + 1 # one additional column since col #0 cannot be modified as the others can
        columns = [f"{col}" for col in range(cols)]
        self.table.configure(columns=columns, displaycolumns=columns)
        

        for col in range(cols):
            if col > 0:
                # set text for the headings
                self.table.heading(columns[col], text=f"{factor_list[col-1]}", anchor='center')
            # set features for columns
            self.table.column(columns[col], width=70, anchor='center', stretch=0)

        # color scheme
        self.table.tag_configure("truth_row", background="white")
        self.table.tag_configure("ordering_row", background="aliceblue")
        
        for row in range(rows-1):
            # fill with zeros
            self.table.insert('', 'end', iid=row, values=['#' + str(row + 1)] + [f"0" for col in range(1,cols)], tag="truth_row")
        # leave bottom row empty
        self.table.insert('', 'end', iid=rows, values=['ordering'] + ["" for col in range(cols)], tag="ordering_row")
        
        self.table.height = rows  # doesn't work??

        # Reconfigure treeview based on the update to resize
        self.table.configure(show='headings')
    
    def create_notebook(self):
        self.notebook = ttk.Notebook(self.root, width=300, height=280)
        self.notebook.grid(row=0, column=0, columnspan=3, sticky="nsew")
        self.root.grid_rowconfigure(0, weight=1)
            
        # Add a table tab to the notebook
        self.create_table_tab(self.notebook)
        
        self.populate_notebook(self.notebook)

    def create_table_tab(self, notebook):
        frame_direct = Frame(notebook)
        self.notebook.add(frame_direct, text="Enter coincidence data")  # Insert at the beginning

        # create label, button and entry-widgets
        label_var = Label(frame_direct, text=f"number of variables")
        label_var.grid(row=0, column=0, pady=10, padx=5)
        
        self.num_variables = StringVar()
        self.num_variables.set(str(default_num_var)) 
        entry_num_var = Entry(frame_direct, width=3, textvariable=self.num_variables)
        entry_num_var.grid(row=0, column=1, sticky=(W, E))
        
        label_conf = Label(frame_direct, text=f"number of configurations")
        label_conf.grid(row=0, column=2, pady=10, padx=5)
        
        self.num_configurations = StringVar()
        self.num_configurations.set(str(default_num_conf))
        entry_num_conf = Entry(frame_direct, width=4, textvariable=self.num_configurations)
        entry_num_conf.grid(row=0, column=3, sticky=(W, E))
        
        button = Button(frame_direct, text="Update", command=self.update)
        button.grid(row=0, column=4, pady=5, padx=5, sticky=W)
        
        label_expl = Label(frame_direct, text=f"Click on a cell to change its value.")
        label_expl.grid(row=2, column=0, columnspan=5, pady=10, padx=5)
        
        # Create a Treeview widget for the table
        self.table = ttk.Treeview(frame_direct, columns=tuple(f"Col{i}" for i in range(1, default_num_var)), show="headings", height = 21)
        
        self.update(new=True)

        # Enable cell editing
        self.table.bind("<ButtonRelease-1>", lambda event: self.on_click(event))
        

        # Pack the table and configure weights
        self.table.grid(row=1, column=0, columnspan=5)

        
    def on_click(self, event):
        region = self.table.identify("region", event.x, event.y)
        if region == "heading":
            def ok(event):
                """Change heading text."""
                self.table.heading(column, text=entry.get())
                if entry in entry_list:
                    entry_list.remove(entry)
                entry.destroy()

            column = self.table.identify_column(event.x) # identify column
            # tree.bbox work sonly with items so we have to get the bbox of the heading differently
            x, y, width, _ = self.table.bbox(self.table.get_children('')[0], column) # get x and width (same as the one of any cell in the column)
            # get vertical coordinates (y1, y2)
            y2 = y
            # get bottom coordinate
            while self.table.identify_region(event.x, y2) != 'heading':  
                y2 -= 1
            # get top coordinate
            y1 = y2
            while self.table.identify_region(event.x, y1) == 'heading':
                y1 -= 1
            height = y2 - y1
            y = y1
            value = self.table.heading(column, 'text')
            
            # display the Entry   
            entry = ttk.Entry(self.table)  # create edition entry
            entry.place(x=x, y=y, width=width, height=height,
                anchor='nw')  # display entry on top of cell
            entry.insert(0, value)  # put former value in entry
            entry.bind('<FocusOut>', lambda e: entry.destroy())  
            entry.bind('<Return>', ok)  # validate with Enter
            entry.focus_set()
            entry_list.append(entry)
            
        elif self.table.selection():
            row = self.table.selection()[0]
            col = self.table.identify_column(event.x)
            
            try:
                index = int(col.split("#")[-1]) - 1
            
                if index > 0:
                    if self.table.item(self.table.focus())['values'][index] == 0:
                        new_value = 1
                    elif self.table.item(self.table.focus())['values'][index] == 1:
                        new_value = 0
                    elif self.table.item(self.table.focus())['values'][index] == "":
                        new_value = "<"
                    elif self.table.item(self.table.focus())['values'][index] == "<":
                        new_value = "<<"
                    elif self.table.item(self.table.focus())['values'][index] == "<<":
                        new_value = ""
                    self.table.set(row, column=col, value=new_value)
            except:
                pass
        

    def populate_notebook(self, notebook):
        # create frames
        frame_csv = Frame(self.notebook, width=300, height=280)
        frame_cna = Frame(self.notebook, width=300, height=280)

        # add frames to notebook
        self.notebook.add(frame_csv, text='Import data from csv')
        self.notebook.add(frame_cna, text='Import data from QCA/CNA')       
        
        # widgets on frame_csv                      
        def select_file():
            filename = filedialog.askopenfilename(title="Select a file", filetypes = (("all files", "*.*"), ("csv files", "*.csv*")))
            if filename:
               self.label_file_csv.configure(text=filename)
               self.label_file_cna.configure(text=filename)
        
        # button to trigger file selection
        select_button = Button(frame_csv, text="Select File", command=select_file)
        select_button.grid(row=0, column=0, pady=5, padx=5, sticky=W)
        
        self.label_file_csv = Label(frame_csv, text="", anchor="w", wraplength=400)
        self.label_file_csv.grid(row=1, column=0, pady=3, padx=5, sticky="nsew")
        
        label_specs = Label(frame_csv, text="Specifications for preparing the csv files:", anchor="w", justify=LEFT)
        label_specs.grid(row=2, column=0, pady=0, padx=5, sticky=W)
        
        desc_csv = " - Enter the variable names into the first row\n - They have to be ordered by ascending constitutive level and if a causal ordering is provided, causally downstream within each level\n - Truth values can be written as \"T\"/\"F\", \"t\"/\"f\" or \"0\"/\"1\"\n - Specify information on the level hierarchy and/or causal ordering in the row below the last row with truth values.\n - Write \"<<\" into the respective column of the first variable belonging to the next higher level\n - Causal ordering can be indicated in the same way by adding \"<\""
        label_desc_csv = Label(frame_csv, text=desc_csv, anchor="w", justify=LEFT, wraplength=400)
        label_desc_csv.grid(row=3, column=0, pady=0, padx=5, columnspan=2, sticky=W)
        
        frame_csv.columnconfigure(1, weight=3)
        
        # to be added -> display data from csv
        # https://www.w3resource.com/python-exercises/tkinter/python-tkinter-dialogs-and-file-handling-exercise-9.php
        
        # widgets on frame_cna
        label_spec_cna = Label(frame_cna, wraplength=400, text="Export the output of QCA or CNA into a text file. Select the file for data import.")
        #label = Label(frame_cna, text=f"{title} - Tab")
        label_spec_cna.grid(row=0, column=0, pady=5, padx=5, columnspan=2, sticky=W)    
        
        # button to trigger file selection
        select_button = Button(frame_cna, text="Select File", command=select_file)
        select_button.grid(row=1, column=0, pady=5, padx=5, sticky=W)    
        
        self.label_file_cna = Label(frame_cna, text="", wraplength=300)
        self.label_file_cna.grid(row=1, column=1, pady=5, padx=5, sticky=W)        
    
    def create_constraints(self):
        frame_constraints = ttk.Frame(self.root)
        frame_constraints.grid(row=1, column=0, columnspan=3, sticky="nsew")
        self.root.grid_rowconfigure(0, weight=1)
            
        # add widgets
        self.populate_constraints(frame_constraints)
    
    def populate_constraints(self, frame):
        # add dummy button
        button_dummy = Button(frame, text="Contraints - to be implemented", state = "disabled")
        button_dummy.grid(row=0, column=2, columnspan=3, pady=10, padx=10, sticky=W)      
    
    def create_run(self):
        frame_run = ttk.Frame(self.root)
        frame_run.grid(row=2, column=0, columnspan=3, sticky="nsew")
        self.root.grid_rowconfigure(0, weight=1)
            
        # add run button
        button_run = Button(frame_run, text="Run mLCA", command=self.run_mlca)
        button_run.grid(row=0, column=0, pady=10, padx=10, sticky=W+E)   
        
        # checkboxes for color-bw mode // simple-complex mode // export list as pdf
        self.var_inus = BooleanVar(value=True)
        checkbox_causal_mode = Checkbutton(frame_run, text='INUS',variable=self.var_inus, onvalue=True, offvalue=False)
        checkbox_causal_mode.grid(row=0, column=1, pady=10, padx=10, sticky=W+E)   

        self.var_color = BooleanVar(value=True)
        checkbox_color = Checkbutton(frame_run, text='color output',variable=self.var_color, onvalue=True, offvalue=False)
        checkbox_color.grid(row=0, column=2, pady=10, padx=10, sticky=W+E)

        self.var_export_list = BooleanVar(value=False)
        checkbox_export_list = Checkbutton(frame_run, text='export list of solutions',variable=self.var_export_list, onvalue=True, offvalue=False)
        checkbox_export_list.grid(row=0, column=3, pady=10, padx=10, sticky=W+E)   

    def export_to_csv(self, file_path):
        if file_path:
            with open(file_path, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
            
                # Write column headers
                headers = [self.table.heading(column)['text'] for column in self.table['columns'][1:]] 
                # do not export the content of first column ([1:]) - it does not contain truth values
                csv_writer.writerow(headers)
            
                # Write data rows
                for row_id in self.table.get_children():
                    row_data = self.table.item(row_id)['values'][1:] # do not export the content of first column ([1:])
                    csv_writer.writerow(row_data)

            print("Data exported to CSV file successfully.")
    
    def run_mlca(self):
        # read checkboxes for mode options
        force_mode = []
        if self.var_inus.get():
            force_mode.append("simple")
        else:
            pass   
        if self.var_color.get():
            force_mode.append("color")
        else:
            force_mode.append("bw")
        if self.var_export_list.get():
            force_mode.append("fulllist")
        else:
            pass
        
        # check which notebook is opened to run with the appropriate data import
        if self.notebook.index("current") == 0:
            # use data from treeview table
            # export it as csv and start mLCA by reading this file
            # get current time for timestamp as file name
            now = datetime.now()
            
            file_path = "dataset_" + now.strftime("%Y_%m_%d_%H_%M_%S") + ".csv" # 
            
            # export as csv
            self.export_to_csv(file_path)
            time.sleep(1) # wait a second for writting csv file to disk
            # run mLCA in csv-mode
            mLCA.main(input_file=file_path, input_type='csv', force_mode=force_mode)
            
        elif self.notebook.index("current") == 1:
            # work with csv file from drive
            mLCA.main(input_file=self.label_file_csv.cget("text"), input_type='csv', force_mode=force_mode)
        elif self.notebook.index("current") == 2:
            # use R-output file
            mLCA.main(input_file=self.label_file_cna.cget("text"), input_type='R', force_mode=force_mode)
        
        # now display pdf
        time.sleep(1) # wait a second
        self.open_file()
        self.display_page()
    
    def create_outer_table(self):
        outer_table_frame = ttk.Frame(self.root)
        outer_table_frame.grid(row=0, column=6, rowspan=3, sticky="nsew")
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(2, weight=1)
        """
        list asf etc. - not yet implemented
        # Create a Treeview widget for the table
        self.outer_table = ttk.Treeview(outer_table_frame, columns=tuple(f"Col{i}" for i in range(1, self.num_columns + 1)), show="headings")

        for i in range(1, self.num_columns + 1):
            column_name = f"Column {i}"
            self.outer_table.heading(f"Col{i}", text=column_name)

        # Populate the table with dummy data
        for i in range(1, self.num_rows + 1):
            values = tuple(f"Data {j}" for j in range(1, self.num_columns + 1))
            self.outer_table.insert("", "end", values=(f"Row {i}", *values))


        # Pack the table and configure weights
        self.outer_table.pack(fill="both", expand=True)
        outer_table_frame.grid_rowconfigure(0, weight=1)
        outer_table_frame.grid_columnconfigure(0, weight=1)
        """
    
    def create_pdf_display(self):
        self.top_frame = ttk.Frame(self.root, width=pdf_width+20, height=pdf_height+20)
        # placing the frame using inside main window using grid()
        self.top_frame.grid(row=0, column=6, pady=10, padx=10)
        # the frame will not propagate
        self.top_frame.grid_propagate(False)
        # creating the bottom frame
        self.bottom_frame = ttk.Frame(self.root, width=pdf_width+20, height=100)
        # placing the frame using inside main window using grid()
        self.bottom_frame.grid(row=1, column=6, pady=10, padx=10)
        # the frame will not propagate
        self.bottom_frame.grid_propagate(False)
        
        # creating a vertical scrollbar
        self.scrolly = Scrollbar(self.top_frame, orient=VERTICAL)
        # adding the scrollbar
        self.scrolly.grid(row=0, column=1, sticky=(N,S))
        # creating a horizontal scrollbar
        self.scrollx = Scrollbar(self.top_frame, orient=HORIZONTAL)
        # adding the scrollbar
        self.scrollx.grid(row=1, column=0, sticky=(W, E))
        
        # creating the canvas for display the PDF pages
        self.output = Canvas(self.top_frame, bg='#FFFFFF', width=560, height=640)
        # inserting both vertical and horizontal scrollbars to the canvas
        self.output.configure(yscrollcommand=self.scrolly.set, xscrollcommand=self.scrollx.set)
        Image.MAX_IMAGE_PIXELS = 10000  # suppress DecompressionBombError for the big image
        # adding the canvas
        self.output.grid(row=0, column=0)
        # configuring the horizontal scrollbar to the canvas
        self.scrolly.configure(command=self.output.yview)
        # configuring the vertical scrollbar to the canvas
        self.scrollx.configure(command=self.output.xview)
        
        # configuring the vertical scrollbar to the canvas
        self.scrollx.configure(command=self.output.xview)
        
        self.upbutton = ttk.Button(self.bottom_frame, text='Previous', command=self.previous_page)
        # adding the button
        self.upbutton.grid(row=0, column=1, padx=(170, 5), pady=8)
        # creating a down button with an icon
        #self.downbutton = ttk.Button(self.bottom_frame, image=self.downarrow, command=self.next_page)
        self.downbutton = ttk.Button(self.bottom_frame, text='Next', command=self.next_page)
        # adding the button
        self.downbutton.grid(row=0, column=2, padx=5, pady=8)
        # creating zoom buttons
        self.zoomoutbutton = ttk.Button(self.bottom_frame, text='Zoom out', command=self.zoom_out)
        # adding the button
        self.zoomoutbutton.grid(row=0, column=3, padx=5, pady=8)
        self.zoominbutton = ttk.Button(self.bottom_frame, text='Zoom in', command=self.zoom_in)
        # adding the button
        self.zoominbutton.grid(row=0, column=4, padx=5, pady=8)
        # setting initial zoom modifier
        self.zoom_modifier = 0
        # label for displaying page numbers
        self.page_label = ttk.Label(self.bottom_frame, text='')
        # adding the label
        self.page_label.grid(row=1, column=2, padx=5, sticky = (W+E))
        
    # function for opening pdf files
    def open_file(self):
        filepath = "output_graph.pdf"
        # checking if the file exists
        if filepath:
            # declaring the path
            self.path = filepath
            # extracting the pdf file from the path
            filename = os.path.basename(self.path)
            # passing the path to PDFMiner 
            self.miner = PDFObject(self.path)
            # getting data and numPages
            numPages = len(self.miner.pdf)

            # setting the current page to 0
            self.current_page = 0
            # checking if numPages exists
            if numPages:
                self.numPages = numPages
                # setting fileopen to True
                self.fileisopen = True
                # calling the display_page() function
                self.display_page()
                # replacing the window title with the PDF document name


    # the function to display the page  
    def display_page(self):
        # checking if numPages is less than current_page and if current_page is less than
        # or equal to 0
        if 0 <= self.current_page < self.numPages:
            # getting the page using get_page() function from miner
            self.img_file = self.miner.get_page(self.current_page,self.zoom_modifier)
            # inserting the page image inside the Canvas
            self.output.create_image(0, 0, anchor='nw', image=self.img_file)
            # the variable to be stringified
            self.stringified_current_page = self.current_page + 1
            # updating the page label with number of pages 
            self.page_label['text'] = str(self.stringified_current_page) + ' of ' + str(self.numPages)
            # creating a region for inserting the page inside the Canvas
            region = self.output.bbox(ALL)
            # making the region to be scrollable
            self.output.configure(scrollregion=region)   
            
    # function for displaying next page
    def next_page(self):
        # checking if file is open
        if self.fileisopen:
            # checking if current_page is less than or equal to numPages-1
            if self.current_page <= self.numPages - 1:
                # updating the page with value 1
                self.current_page += 1
                # displaying the new page
                aux = self.zoom_modifier
                self.zoom_modifier = 0
                self.display_page() 
                self.zoom_modifier = aux 

    # function for displaying the previous page        
    def previous_page(self):
        # checking if fileisopen
        if self.fileisopen:
            # checking if current_page is greater than 0
            if self.current_page > 0:
                # decrementing the current_page by 1
                self.current_page -= 1
                # displaying the previous page
                aux = self.zoom_modifier
                self.zoom_modifier = 0
                self.display_page() 
                self.zoom_modifier = aux                           
    
    def zoom_in(self):
        try:
            if self.fileisopen:
                if self.zoom_modifier == 0:
                    self.zoom_modifier = 1
                    self.display_page()     
                elif self.zoom_modifier == -1:
                    self.zoom_modifier = 0
                    self.display_page()
        except:
            pass
                
    
    def zoom_out(self):
        try:
            if self.fileisopen:
                if self.zoom_modifier == 1:
                    self.zoom_modifier = 0
                    self.display_page()
                elif self.zoom_modifier == 0:
                    self.zoom_modifier = -1
                    self.display_page()
        except:
            pass
        
        
    
if __name__ == "__main__":
    root = Tk()
    app = NotebookGridApp(root)
    root.mainloop()
