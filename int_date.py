from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QComboBox, QPushButton
from PyQt5.QtCore import Qt
import sys
import os
from functools import partial
from logistics import *

class App(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = 'Logistics App'
        self.left = 10
        self.top = 10
        self.width = 400
        self.height = 300
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Add a label
        label = QLabel('Select a date:', self)
        label.setAlignment(Qt.AlignCenter)
        label.setGeometry(0, 20, self.width, 30)

        # Add a combobox
        self.combobox = QComboBox(self)
        self.combobox.setGeometry(50, 80, 300, 30)

        # Add a button
        button = QPushButton('Get Results', self)
        button.setGeometry(150, 150, 100, 30)
        button.clicked.connect(self.getResults)

    def fillComboBox(self):
        # Get all dates from the Outputs directory
        options = os.listdir(os.path.join("Outputs", ))
        self.combobox.addItems(options)

    def getResults(self):
        # Get the selected date from the combobox
        date = self.combobox.currentText()

        # Call wanted_date function with the selected date
        directory_read = os.path.join("Outputs", date)
        result = wanted_date(directory_read)

        # Display the result in a label
        label = QLabel(result, self)
        label.setGeometry(0, 200, self.width, 50)
        label.setAlignment(Qt.AlignCenter)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    ex.fillComboBox()
    ex.show()
    sys.exit(app.exec_())
