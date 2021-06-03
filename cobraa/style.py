from ROOT import TStyle,gROOT

def setStyle():

    cobraaStyle = TStyle("cobraa","Style for Cobraa plots")
    tsize = 0.05
    toffset = 1.4

    cobraaStyle.SetPadTopMargin(0.1)
    cobraaStyle.SetPadRightMargin(0.1)
    cobraaStyle.SetPadBottomMargin(0.16)
    cobraaStyle.SetPadLeftMargin(0.16)


    cobraaStyle.SetTextSize(tsize)
    cobraaStyle.SetLabelSize(tsize,'x')
    cobraaStyle.SetTitleSize(tsize,'x')
    cobraaStyle.SetLabelSize(tsize,'y')
    cobraaStyle.SetTitleSize(tsize,'y')
    cobraaStyle.SetLabelSize(tsize,'z')
    cobraaStyle.SetTitleSize(tsize,'z')
    cobraaStyle.SetTitleXOffset(toffset)
    cobraaStyle.SetTitleYOffset(toffset)

    cobraaStyle.SetMarkerStyle(20)
    cobraaStyle.SetMarkerSize(1.2)
    cobraaStyle.SetHistLineWidth(2)
    cobraaStyle.SetLineStyleString(2,'[12 12]')


    cobraaStyle.SetPadTopMargin(0.05)
    cobraaStyle.SetPadRightMargin(0.05)
    cobraaStyle.SetPadBottomMargin(0.16)
    cobraaStyle.SetPadLeftMargin(0.16)

    gROOT.SetStyle("cobraa")

    return cobraaStyle
