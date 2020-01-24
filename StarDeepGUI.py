# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:10:15 2020

@author: Red
"""

import sys
from PyQt5 import QtWidgets, uic
  
app = QtWidgets.QApplication(sys.argv)

window = uic.loadUi("StarDeepUI.ui")
window.show()
app.exec()