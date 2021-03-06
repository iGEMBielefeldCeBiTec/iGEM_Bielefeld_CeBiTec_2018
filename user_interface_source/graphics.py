# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'siRNA_Graphics.ui'
#
# Created by: PyQt4 UI code generator 4.12.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(978, 847)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8("siRCon_Logo_1.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setAutoFillBackground(False)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setToolTip(_fromUtf8(""))
        self.tabWidget.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.tabWidget.setTabShape(QtGui.QTabWidget.Rounded)
        self.tabWidget.setElideMode(QtCore.Qt.ElideNone)
        self.tabWidget.setDocumentMode(False)
        self.tabWidget.setTabsClosable(False)
        self.tabWidget.setMovable(False)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab = QtGui.QWidget()
        self.tab.setEnabled(True)
        self.tab.setObjectName(_fromUtf8("tab"))
        self.verticalLayout_4 = QtGui.QVBoxLayout(self.tab)
        self.verticalLayout_4.setMargin(0)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setContentsMargins(10, 20, 10, -1)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label = QtGui.QLabel(self.tab)
        self.label.setToolTip(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout_3.addWidget(self.label)
        self.geneSeq1 = QtGui.QTextEdit(self.tab)
        self.geneSeq1.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.geneSeq1.setObjectName(_fromUtf8("geneSeq1"))
        self.verticalLayout_3.addWidget(self.geneSeq1)
        self.verticalLayout_4.addLayout(self.verticalLayout_3)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setContentsMargins(10, 20, 10, 20)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.gridLayout_4 = QtGui.QGridLayout()
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.checkBox = QtGui.QCheckBox(self.tab)
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.gridLayout_4.addWidget(self.checkBox, 0, 0, 1, 1)
        self.horizontalLayout_2.addLayout(self.gridLayout_4)
        self.gridLayout_6 = QtGui.QGridLayout()
        self.gridLayout_6.setObjectName(_fromUtf8("gridLayout_6"))
        self.label_8 = QtGui.QLabel(self.tab)
        self.label_8.setText(_fromUtf8(""))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_6.addWidget(self.label_8, 0, 0, 1, 1)
        self.horizontalLayout_2.addLayout(self.gridLayout_6)
        self.gridLayout_5 = QtGui.QGridLayout()
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.getRNA1 = QtGui.QPushButton(self.tab)
        self.getRNA1.setEnabled(True)
        self.getRNA1.setMaximumSize(QtCore.QSize(174, 27))
        self.getRNA1.setSizeIncrement(QtCore.QSize(0, 0))
        self.getRNA1.setToolTip(_fromUtf8(""))
        self.getRNA1.setObjectName(_fromUtf8("getRNA1"))
        self.gridLayout_5.addWidget(self.getRNA1, 0, 0, 1, 1)
        self.horizontalLayout_2.addLayout(self.gridLayout_5)
        self.verticalLayout_4.addLayout(self.horizontalLayout_2)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setContentsMargins(10, -1, 10, -1)
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.label_3 = QtGui.QLabel(self.tab)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_5.addWidget(self.label_3)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout_16 = QtGui.QVBoxLayout()
        self.verticalLayout_16.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.verticalLayout_16.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout_16.setSpacing(6)
        self.verticalLayout_16.setObjectName(_fromUtf8("verticalLayout_16"))
        self.label_15 = QtGui.QLabel(self.tab)
        self.label_15.setMaximumSize(QtCore.QSize(99, 16777215))
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.verticalLayout_16.addWidget(self.label_15)
        self.output_text1_3 = QtGui.QTextBrowser(self.tab)
        self.output_text1_3.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.output_text1_3.setObjectName(_fromUtf8("output_text1_3"))
        self.verticalLayout_16.addWidget(self.output_text1_3)
        self.gridLayout.addLayout(self.verticalLayout_16, 0, 1, 1, 1)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        self.verticalLayout_14 = QtGui.QVBoxLayout()
        self.verticalLayout_14.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout_14.setObjectName(_fromUtf8("verticalLayout_14"))
        self.label_2 = QtGui.QLabel(self.tab)
        self.label_2.setMaximumSize(QtCore.QSize(279, 16777215))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_14.addWidget(self.label_2)
        self.output_text1_1 = QtGui.QTextBrowser(self.tab)
        self.output_text1_1.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.output_text1_1.setObjectName(_fromUtf8("output_text1_1"))
        self.verticalLayout_14.addWidget(self.output_text1_1)
        self.horizontalLayout_7.addLayout(self.verticalLayout_14)
        self.verticalLayout_15 = QtGui.QVBoxLayout()
        self.verticalLayout_15.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout_15.setObjectName(_fromUtf8("verticalLayout_15"))
        self.label_11 = QtGui.QLabel(self.tab)
        self.label_11.setMaximumSize(QtCore.QSize(279, 16777215))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.verticalLayout_15.addWidget(self.label_11)
        self.output_text1_2 = QtGui.QTextBrowser(self.tab)
        self.output_text1_2.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.output_text1_2.setObjectName(_fromUtf8("output_text1_2"))
        self.verticalLayout_15.addWidget(self.output_text1_2)
        self.horizontalLayout_7.addLayout(self.verticalLayout_15)
        self.gridLayout.addLayout(self.horizontalLayout_7, 0, 0, 1, 1)
        self.gridLayout.setColumnStretch(0, 9)
        self.gridLayout.setColumnStretch(1, 2)
        self.verticalLayout_5.addLayout(self.gridLayout)
        self.verticalLayout_4.addLayout(self.verticalLayout_5)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setMargin(10)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.gridLayout_7 = QtGui.QGridLayout()
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        self.label_9 = QtGui.QLabel(self.tab)
        self.label_9.setText(_fromUtf8(""))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_7.addWidget(self.label_9, 0, 1, 1, 1)
        self.horizontalLayout_3.addLayout(self.gridLayout_7)
        self.gridLayout_9 = QtGui.QGridLayout()
        self.gridLayout_9.setObjectName(_fromUtf8("gridLayout_9"))
        self.scaffold_1 = QtGui.QCheckBox(self.tab)
        self.scaffold_1.setObjectName(_fromUtf8("scaffold_1"))
        self.gridLayout_9.addWidget(self.scaffold_1, 0, 0, 1, 1)
        self.horizontalLayout_3.addLayout(self.gridLayout_9)
        self.gridLayout_10 = QtGui.QGridLayout()
        self.gridLayout_10.setObjectName(_fromUtf8("gridLayout_10"))
        self.save_button1 = QtGui.QPushButton(self.tab)
        self.save_button1.setEnabled(True)
        self.save_button1.setMaximumSize(QtCore.QSize(235, 27))
        self.save_button1.setSizeIncrement(QtCore.QSize(0, 0))
        self.save_button1.setObjectName(_fromUtf8("save_button1"))
        self.gridLayout_10.addWidget(self.save_button1, 0, 0, 1, 1)
        self.horizontalLayout_3.addLayout(self.gridLayout_10)
        self.verticalLayout_4.addLayout(self.horizontalLayout_3)
        self.textBrowser_2 = QtGui.QTextBrowser(self.tab)
        self.textBrowser_2.setObjectName(_fromUtf8("textBrowser_2"))
        self.verticalLayout_4.addWidget(self.textBrowser_2)
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.verticalLayout_9 = QtGui.QVBoxLayout(self.tab_2)
        self.verticalLayout_9.setMargin(0)
        self.verticalLayout_9.setObjectName(_fromUtf8("verticalLayout_9"))
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setContentsMargins(10, 20, 10, -1)
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.label_12 = QtGui.QLabel(self.tab_2)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.verticalLayout_6.addWidget(self.label_12)
        self.geneSeq2 = QtGui.QTextEdit(self.tab_2)
        self.geneSeq2.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.geneSeq2.setObjectName(_fromUtf8("geneSeq2"))
        self.verticalLayout_6.addWidget(self.geneSeq2)
        self.verticalLayout_9.addLayout(self.verticalLayout_6)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setContentsMargins(10, 20, 10, 20)
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.gridLayout_11 = QtGui.QGridLayout()
        self.gridLayout_11.setContentsMargins(-1, -1, 1, -1)
        self.gridLayout_11.setObjectName(_fromUtf8("gridLayout_11"))
        self.checkBox2 = QtGui.QCheckBox(self.tab_2)
        self.checkBox2.setObjectName(_fromUtf8("checkBox2"))
        self.gridLayout_11.addWidget(self.checkBox2, 0, 0, 1, 1)
        self.horizontalLayout_4.addLayout(self.gridLayout_11)
        self.gridLayout_15 = QtGui.QGridLayout()
        self.gridLayout_15.setObjectName(_fromUtf8("gridLayout_15"))
        self.label_4 = QtGui.QLabel(self.tab_2)
        self.label_4.setText(_fromUtf8(""))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_15.addWidget(self.label_4, 0, 0, 1, 1)
        self.horizontalLayout_4.addLayout(self.gridLayout_15)
        self.gridLayout_16 = QtGui.QGridLayout()
        self.gridLayout_16.setObjectName(_fromUtf8("gridLayout_16"))
        self.getRNA2 = QtGui.QPushButton(self.tab_2)
        self.getRNA2.setMaximumSize(QtCore.QSize(174, 27))
        self.getRNA2.setObjectName(_fromUtf8("getRNA2"))
        self.gridLayout_16.addWidget(self.getRNA2, 0, 0, 1, 1)
        self.horizontalLayout_4.addLayout(self.gridLayout_16)
        self.verticalLayout_9.addLayout(self.horizontalLayout_4)
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setContentsMargins(10, -1, 10, -1)
        self.verticalLayout_8.setObjectName(_fromUtf8("verticalLayout_8"))
        self.label_5 = QtGui.QLabel(self.tab_2)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.verticalLayout_8.addWidget(self.label_5)
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName(_fromUtf8("horizontalLayout_8"))
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName(_fromUtf8("horizontalLayout_9"))
        self.verticalLayout_17 = QtGui.QVBoxLayout()
        self.verticalLayout_17.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout_17.setObjectName(_fromUtf8("verticalLayout_17"))
        self.label_17 = QtGui.QLabel(self.tab_2)
        self.label_17.setMaximumSize(QtCore.QSize(279, 16777215))
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.verticalLayout_17.addWidget(self.label_17)
        self.output_text2_1 = QtGui.QTextBrowser(self.tab_2)
        self.output_text2_1.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.output_text2_1.setObjectName(_fromUtf8("output_text2_1"))
        self.verticalLayout_17.addWidget(self.output_text2_1)
        self.horizontalLayout_9.addLayout(self.verticalLayout_17)
        self.verticalLayout_19 = QtGui.QVBoxLayout()
        self.verticalLayout_19.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout_19.setObjectName(_fromUtf8("verticalLayout_19"))
        self.label_18 = QtGui.QLabel(self.tab_2)
        self.label_18.setMaximumSize(QtCore.QSize(279, 16777215))
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.verticalLayout_19.addWidget(self.label_18)
        self.output_text2_2 = QtGui.QTextBrowser(self.tab_2)
        self.output_text2_2.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.output_text2_2.setObjectName(_fromUtf8("output_text2_2"))
        self.verticalLayout_19.addWidget(self.output_text2_2)
        self.horizontalLayout_9.addLayout(self.verticalLayout_19)
        self.horizontalLayout_8.addLayout(self.horizontalLayout_9)
        self.verticalLayout_23 = QtGui.QVBoxLayout()
        self.verticalLayout_23.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        self.verticalLayout_23.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout_23.setSpacing(6)
        self.verticalLayout_23.setObjectName(_fromUtf8("verticalLayout_23"))
        self.label_19 = QtGui.QLabel(self.tab_2)
        self.label_19.setMaximumSize(QtCore.QSize(99, 16777215))
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.verticalLayout_23.addWidget(self.label_19)
        self.output_text2_3 = QtGui.QTextBrowser(self.tab_2)
        self.output_text2_3.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.output_text2_3.setObjectName(_fromUtf8("output_text2_3"))
        self.verticalLayout_23.addWidget(self.output_text2_3)
        self.horizontalLayout_8.addLayout(self.verticalLayout_23)
        self.horizontalLayout_8.setStretch(0, 9)
        self.horizontalLayout_8.setStretch(1, 2)
        self.verticalLayout_8.addLayout(self.horizontalLayout_8)
        self.verticalLayout_9.addLayout(self.verticalLayout_8)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setMargin(10)
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.gridLayout_17 = QtGui.QGridLayout()
        self.gridLayout_17.setObjectName(_fromUtf8("gridLayout_17"))
        self.scaffold_3 = QtGui.QCheckBox(self.tab_2)
        self.scaffold_3.setObjectName(_fromUtf8("scaffold_3"))
        self.gridLayout_17.addWidget(self.scaffold_3, 0, 0, 1, 1)
        self.horizontalLayout_5.addLayout(self.gridLayout_17)
        self.gridLayout_18 = QtGui.QGridLayout()
        self.gridLayout_18.setObjectName(_fromUtf8("gridLayout_18"))
        self.scaffold_2 = QtGui.QCheckBox(self.tab_2)
        self.scaffold_2.setObjectName(_fromUtf8("scaffold_2"))
        self.gridLayout_18.addWidget(self.scaffold_2, 0, 0, 1, 1)
        self.horizontalLayout_5.addLayout(self.gridLayout_18)
        self.gridLayout_20 = QtGui.QGridLayout()
        self.gridLayout_20.setObjectName(_fromUtf8("gridLayout_20"))
        self.save_button2 = QtGui.QPushButton(self.tab_2)
        self.save_button2.setEnabled(True)
        self.save_button2.setMaximumSize(QtCore.QSize(235, 27))
        self.save_button2.setSizeIncrement(QtCore.QSize(0, 0))
        self.save_button2.setObjectName(_fromUtf8("save_button2"))
        self.gridLayout_20.addWidget(self.save_button2, 0, 0, 1, 1)
        self.horizontalLayout_5.addLayout(self.gridLayout_20)
        self.verticalLayout_9.addLayout(self.horizontalLayout_5)
        self.textBrowser_3 = QtGui.QTextBrowser(self.tab_2)
        self.textBrowser_3.setObjectName(_fromUtf8("textBrowser_3"))
        self.verticalLayout_9.addWidget(self.textBrowser_3)
        self.tabWidget.addTab(self.tab_2, _fromUtf8(""))
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.verticalLayout_11 = QtGui.QVBoxLayout(self.tab_3)
        self.verticalLayout_11.setMargin(0)
        self.verticalLayout_11.setObjectName(_fromUtf8("verticalLayout_11"))
        self.verticalLayout_18 = QtGui.QVBoxLayout()
        self.verticalLayout_18.setObjectName(_fromUtf8("verticalLayout_18"))
        self.verticalLayout_20 = QtGui.QVBoxLayout()
        self.verticalLayout_20.setContentsMargins(10, 10, 10, -1)
        self.verticalLayout_20.setObjectName(_fromUtf8("verticalLayout_20"))
        self.label_20 = QtGui.QLabel(self.tab_3)
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.verticalLayout_20.addWidget(self.label_20)
        self.geneSeq3 = QtGui.QTextEdit(self.tab_3)
        self.geneSeq3.setInputMethodHints(QtCore.Qt.ImhNone)
        self.geneSeq3.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.geneSeq3.setTabChangesFocus(False)
        self.geneSeq3.setObjectName(_fromUtf8("geneSeq3"))
        self.verticalLayout_20.addWidget(self.geneSeq3)
        self.verticalLayout_21 = QtGui.QVBoxLayout()
        self.verticalLayout_21.setObjectName(_fromUtf8("verticalLayout_21"))
        self.label_21 = QtGui.QLabel(self.tab_3)
        self.label_21.setToolTip(_fromUtf8(""))
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.verticalLayout_21.addWidget(self.label_21)
        self.siRNASeq = QtGui.QTextEdit(self.tab_3)
        self.siRNASeq.setMaximumSize(QtCore.QSize(16777215, 31))
        self.siRNASeq.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.siRNASeq.setObjectName(_fromUtf8("siRNASeq"))
        self.verticalLayout_21.addWidget(self.siRNASeq)
        self.verticalLayout_20.addLayout(self.verticalLayout_21)
        self.verticalLayout_18.addLayout(self.verticalLayout_20)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setContentsMargins(-1, 20, -1, 20)
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.gridLayout_37 = QtGui.QGridLayout()
        self.gridLayout_37.setContentsMargins(10, -1, -1, -1)
        self.gridLayout_37.setObjectName(_fromUtf8("gridLayout_37"))
        self.method = QtGui.QComboBox(self.tab_3)
        self.method.setMaximumSize(QtCore.QSize(181, 27))
        self.method.setObjectName(_fromUtf8("method"))
        self.method.addItem(_fromUtf8(""))
        self.method.addItem(_fromUtf8(""))
        self.method.addItem(_fromUtf8(""))
        self.gridLayout_37.addWidget(self.method, 0, 0, 1, 1)
        self.horizontalLayout_10.addLayout(self.gridLayout_37)
        self.gridLayout_39 = QtGui.QGridLayout()
        self.gridLayout_39.setContentsMargins(10, -1, -1, -1)
        self.gridLayout_39.setObjectName(_fromUtf8("gridLayout_39"))
        self.checkBox3 = QtGui.QCheckBox(self.tab_3)
        self.checkBox3.setObjectName(_fromUtf8("checkBox3"))
        self.gridLayout_39.addWidget(self.checkBox3, 0, 0, 1, 1)
        self.horizontalLayout_10.addLayout(self.gridLayout_39)
        self.gridLayout_40 = QtGui.QGridLayout()
        self.gridLayout_40.setObjectName(_fromUtf8("gridLayout_40"))
        self.checksiRNA = QtGui.QPushButton(self.tab_3)
        self.checksiRNA.setEnabled(True)
        self.checksiRNA.setMaximumSize(QtCore.QSize(233, 27))
        self.checksiRNA.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
        self.checksiRNA.setCheckable(False)
        self.checksiRNA.setObjectName(_fromUtf8("checksiRNA"))
        self.gridLayout_40.addWidget(self.checksiRNA, 1, 0, 1, 1)
        self.horizontalLayout_10.addLayout(self.gridLayout_40)
        self.verticalLayout_18.addLayout(self.horizontalLayout_10)
        self.verticalLayout_22 = QtGui.QVBoxLayout()
        self.verticalLayout_22.setContentsMargins(10, -1, 10, -1)
        self.verticalLayout_22.setObjectName(_fromUtf8("verticalLayout_22"))
        self.label_23 = QtGui.QLabel(self.tab_3)
        self.label_23.setObjectName(_fromUtf8("label_23"))
        self.verticalLayout_22.addWidget(self.label_23)
        self.output_text3 = QtGui.QTextBrowser(self.tab_3)
        self.output_text3.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.output_text3.setObjectName(_fromUtf8("output_text3"))
        self.verticalLayout_22.addWidget(self.output_text3)
        self.verticalLayout_18.addLayout(self.verticalLayout_22)
        self.horizontalLayout_11 = QtGui.QHBoxLayout()
        self.horizontalLayout_11.setContentsMargins(10, 10, 10, -1)
        self.horizontalLayout_11.setObjectName(_fromUtf8("horizontalLayout_11"))
        self.gridLayout_41 = QtGui.QGridLayout()
        self.gridLayout_41.setObjectName(_fromUtf8("gridLayout_41"))
        self.label_13 = QtGui.QLabel(self.tab_3)
        self.label_13.setText(_fromUtf8(""))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_41.addWidget(self.label_13, 0, 0, 1, 1)
        self.horizontalLayout_11.addLayout(self.gridLayout_41)
        self.gridLayout_42 = QtGui.QGridLayout()
        self.gridLayout_42.setObjectName(_fromUtf8("gridLayout_42"))
        self.gridLayout_43 = QtGui.QGridLayout()
        self.gridLayout_43.setObjectName(_fromUtf8("gridLayout_43"))
        self.label_14 = QtGui.QLabel(self.tab_3)
        self.label_14.setText(_fromUtf8(""))
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_43.addWidget(self.label_14, 0, 0, 1, 1)
        self.gridLayout_42.addLayout(self.gridLayout_43, 0, 0, 1, 1)
        self.horizontalLayout_11.addLayout(self.gridLayout_42)
        self.gridLayout_44 = QtGui.QGridLayout()
        self.gridLayout_44.setContentsMargins(-1, -1, -1, 10)
        self.gridLayout_44.setObjectName(_fromUtf8("gridLayout_44"))
        self.save_button3 = QtGui.QPushButton(self.tab_3)
        self.save_button3.setMaximumSize(QtCore.QSize(235, 27))
        self.save_button3.setObjectName(_fromUtf8("save_button3"))
        self.gridLayout_44.addWidget(self.save_button3, 0, 0, 1, 1)
        self.horizontalLayout_11.addLayout(self.gridLayout_44)
        self.verticalLayout_18.addLayout(self.horizontalLayout_11)
        self.verticalLayout_11.addLayout(self.verticalLayout_18)
        self.tabWidget.addTab(self.tab_3, _fromUtf8(""))
        self.tab_4 = QtGui.QWidget()
        self.tab_4.setObjectName(_fromUtf8("tab_4"))
        self.verticalLayout_12 = QtGui.QVBoxLayout(self.tab_4)
        self.verticalLayout_12.setMargin(0)
        self.verticalLayout_12.setObjectName(_fromUtf8("verticalLayout_12"))
        self.verticalLayout_7 = QtGui.QVBoxLayout()
        self.verticalLayout_7.setObjectName(_fromUtf8("verticalLayout_7"))
        self.textBrowser = QtGui.QTextBrowser(self.tab_4)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.verticalLayout_7.addWidget(self.textBrowser)
        self.verticalLayout_12.addLayout(self.verticalLayout_7)
        self.tabWidget.addTab(self.tab_4, _fromUtf8(""))
        self.verticalLayout.addWidget(self.tabWidget)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setEnabled(False)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 978, 25))
        self.menubar.setNativeMenuBar(True)
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "siRCon", None))
        self.label.setText(_translate("MainWindow", "Gene sequence 5\' -> 3\':", None))
        self.geneSeq1.setProperty("placeholderText", _translate("MainWindow", "Enter gene sequence", None))
        self.checkBox.setText(_translate("MainWindow", "Tace vector system", None))
        self.getRNA1.setText(_translate("MainWindow", "Get siRNA", None))
        self.label_3.setText(_translate("MainWindow", "Results:", None))
        self.label_15.setText(_translate("MainWindow", "Probability", None))
        self.label_2.setText(_translate("MainWindow", "Sense siRNA 5\' -> 3\'", None))
        self.label_11.setText(_translate("MainWindow", "Antisense siRNA 5\' -> 3\'", None))
        self.scaffold_1.setText(_translate("MainWindow", "Save with MicC scaffold", None))
        self.save_button1.setText(_translate("MainWindow", "Save results", None))
        self.textBrowser_2.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt; font-weight:600;\">Tipp (Only when not using the Tace vector system)</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt; font-weight:600;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">In order for the siRNA to bind more strongly to the target sequence, an Hfq binding sequence (here the MicC Hfq sequence) should be added at the 3\' end:</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">sense siRNA + TTTCTGTTGGGCCATTGCATTGCCACTGATTTTCCAACATATAAAAAGACAAGCCCGAACAGTCGTCCGGGCTTTTTTTCTCGAG</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">CTCGAGAAAAAAAGCCCGGACGACTGTTCGGGCTTGTCTTTTTATATGTTGGAAAATCAGTGGCAATGCAATGGCCCAACAGAAA + antisense siRNA</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "siRNA for RNAi", None))
        self.label_12.setToolTip(_translate("MainWindow", "Mark nucleotides with the cursor", None))
        self.label_12.setText(_translate("MainWindow", "Gene sequence 5\' -> 3\':", None))
        self.geneSeq2.setProperty("placeholderText", _translate("MainWindow", "Enter gene sequence", None))
        self.checkBox2.setText(_translate("MainWindow", "Tace vector system", None))
        self.getRNA2.setText(_translate("MainWindow", "Get siRNA", None))
        self.label_5.setText(_translate("MainWindow", "Results:", None))
        self.label_17.setText(_translate("MainWindow", "Sense siRNA 5\' -> 3\'", None))
        self.label_18.setText(_translate("MainWindow", "Antisense siRNA 5\' -> 3\'", None))
        self.label_19.setText(_translate("MainWindow", "Probability", None))
        self.scaffold_3.setText(_translate("MainWindow", "Save with MicC scaffold", None))
        self.scaffold_2.setText(_translate("MainWindow", "Save with OmpA scaffold", None))
        self.save_button2.setText(_translate("MainWindow", "Save results", None))
        self.textBrowser_3.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt; font-weight:600;\">Tipp (Only when not using the Tace vector system)</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt; font-weight:600;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">To increase the stability of the siRNA, the OmpA 5\' UTR can be added to the 5\' end of the siRNA:</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">GCCAGGGGTGCTCGGCATAAGCCGAAGATATCGGTAGAGTTAATATTGAGCAGATCCCCCGGTGAAGGATTTAACCGTGTTATCTCGTTGGAGATATTCATGGCGTATTTTGGATGA + sense siRNA </span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">antisense siRNA + TCATCCAAAATACGCCATGAATATCTCCAACGAGATAACACGGTTAAATCCTTCACCGGGGGATCTGCTCAATATTAACTCTACCGATATCTTCGGCTTATGCCGAGCACCCCTGGC</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">In order for the siRNA to bind more strongly to the target sequence, an Hfq binding sequence (here the MicC Hfq sequence) should be added at the 3\' end:</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">sense siRNA + TTTCTGTTGGGCCATTGCATTGCCACTGATTTTCCAACATATAAAAAGACAAGCCCGAACAGTCGTCCGGGCTTTTTTTCTCGAG</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">CTCGAGAAAAAAAGCCCGGACGACTGTTCGGGCTTGTCTTTTTATATGTTGGAAAATCAGTGGCAATGCAATGGCCCAACAGAAA + antisense siRNA</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "siRNA", None))
        self.label_20.setText(_translate("MainWindow", "Gene sequence 5\' -> 3\':", None))
        self.geneSeq3.setProperty("placeholderText", _translate("MainWindow", "Enter gene sequence", None))
        self.label_21.setText(_translate("MainWindow", "siRNA sequence 5\' -> 3\'", None))
        self.siRNASeq.setProperty("placeholderText", _translate("MainWindow", "Enter siRNA sequence", None))
        self.method.setItemText(0, _translate("MainWindow", "Choose design method", None))
        self.method.setItemText(1, _translate("MainWindow", "RNAi", None))
        self.method.setItemText(2, _translate("MainWindow", "siRNA", None))
        self.checkBox3.setText(_translate("MainWindow", "Tace vector system", None))
        self.checksiRNA.setText(_translate("MainWindow", "Check siRNA", None))
        self.label_23.setText(_translate("MainWindow", "Results:", None))
        self.save_button3.setText(_translate("MainWindow", "Save results", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Check siRNA", None))
        self.textBrowser.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:16pt; font-weight:600;\">Usage</span></p>\n"
"<p align=\"center\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-size:12pt; font-weight:600;\">siRNA for RNAi Module</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 1</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Insert gene sequence</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 2</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Decide if Tace vector system is used (optionally)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 3</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Construction of siRNAs</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 4</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> View resulting siRNAs and their corresponding probability</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 5</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> If Tace system is not used, decide whether MicC scaffold should beadded to the siRNA sequences (optionally)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 6 </span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\">Save results as FASTA file</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'DejaVu Sans, sans-serif\'; font-size:12pt; font-weight:600;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-size:12pt; font-weight:600;\">siRNA Module</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 1</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Insert gene sequence</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 2</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Decide if Tace vector system is used (optionally)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 3</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Construction of siRNAs</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 4</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> View resulting siRNAs and their corresponding probability</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 5</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> If Tace system is not used, decide whether MicC scaffold should be added to the siRNA sequences (optionally)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 6</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> If Tace system is not used, decide whether OmpA scaffold should be added to the siRNA sequences (optionally)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 7</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Save results as FASTA file</span></p>\n"
"<p align=\"justify\" style=\"-qt-paragraph-type:empty; margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'DejaVu Sans, sans-serif\';\"><br /></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">ckeck siRNA (Alpha version)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 1</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Insert gene sequences</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 2</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Insert siRNA sequence</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 3</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Choose the method the siRNA should be used for (siRNA for RNAi or siRNA for silencing)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 4</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Check if the siRNA is used in our Tace system (optionally)</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 5</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Validation of entered siRNA for given target gene sequences</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 6</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> View Results</span></p>\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'DejaVu Sans, sans-serif\'; font-weight:600;\">Step 7</span><span style=\" font-family:\'DejaVu Sans, sans-serif\';\"> Save Results (optionally)</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p align=\"center\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;\"><br /></p>\n"
"<p align=\"center\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;\"><br /></p>\n"
"<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:16pt; font-weight:600;\">Copyright and Licensing</span></p>\n"
"<p align=\"center\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:16pt; font-weight:600;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\"> siRCon - A siRNA Constructor</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Copyright (C) 2018  Vanessa Krämer</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    This program is free software: you can redistribute it and/or modify</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    it under the terms of the GNU General Public License as published by</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    the Free Software Foundation, either version 3 of the License, or</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    (at your option) any later version.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    This program is distributed in the hope that it will be useful,</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    but WITHOUT ANY WARRANTY; without even the implied warranty of</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    GNU General Public License for more details.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:12pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    You should have received a copy of the GNU General Public License</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">    along with this program.  If not, see &lt;https://www.gnu.org/licenses/&gt;.</span></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), _translate("MainWindow", "About", None))

