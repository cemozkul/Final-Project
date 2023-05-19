import sys
import tkinter as tk
    
def help():
    # Open the file and read its contents
    with open("business_definition.txt", "r") as f:
        help_text = f.read()
    
    # Create a new window to show the help text
    help_window = tk.Toplevel()
    help_window.title("Help")
    help_window.geometry("500x400")

    # Add a label widget to show the help text
    help_label = tk.Label(help_window, text=help_text, wraplength=400)
    help_label.pack(expand=True, fill=tk.BOTH)

    # Add a button to close the window
    close_button = tk.Button(help_window, text="Close", command=help_window.destroy)
    close_button.pack(pady=5)

def inputs():
    # Open the file and read its contents
    with open("inputs.txt", "r") as f:
        inputs_text = f.read()
    
    # Create a new window to show the help text
    help_window = tk.Toplevel()
    help_window.title("Help")
    help_window.geometry("500x400")

    # Add a label widget to show the help text
    help_label = tk.Label(help_window, text=inputs_text, wraplength=400)
    help_label.pack(expand=True, fill=tk.BOTH)

    # Add a button to close the window
    close_button = tk.Button(help_window, text="Close", command=help_window.destroy)
    close_button.pack(pady=5)
    
def outputs():
    # Open the file and read its contents
    with open("outputs.txt", "r") as f:
        outputs_text = f.read()
    
    # Create a new window to show the help text
    help_window = tk.Toplevel()
    help_window.title("Help")
    help_window.geometry("500x400")

    # Add a label widget to show the help text
    help_label = tk.Label(help_window, text=outputs_text, wraplength=400)
    help_label.pack(expand=True, fill=tk.BOTH)

    # Add a button to close the window
    close_button = tk.Button(help_window, text="Close", command=help_window.destroy)
    close_button.pack(pady=5)
    
      

