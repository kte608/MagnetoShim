import math

def writepsA3(fileout,Z,R,numZ,numPhi,GrayScale):
    """
    Create and write the postscript code that print the InkShimMatrix
    in three A3 pages.
    Name of the .ps file is what had been put instead of fileout.
    """

    ### size of boxes ###
    dx=Z/float(numZ-1)             #meters 
    dy=2.0*math.pi*R/float(numPhi) #meters 
    dx=int(dx*72/0.0254)           #points = 1/72 inch
    dy=int(dy*72/0.0254)           #points

    ### write the ps file ###
    fh=open(fileout,"w+")
    fh.write("%Postscript to print InkShimMatrix\n")
    
    fh.write("%----Variables----\n")
    fh.write("/dx "+dx.__str__()+" def\n")
    fh.write("/dy "+dy.__str__()+" def\n")
    fh.write("/inch {72 mul}def\n")
    fh.write("/pagewidth 842 def\n")        #A3 format 842*1190 points
    fh.write("/pageheight 1166 def \n") 
    fh.write("/rows 3 def \n")              #number of pages 
    fh.write("/columns 1 def \n")
    fh.write("\n")
        
    fh.write("%----Procedures----\n")
    fh.write("\n")
    fh.write("/TheBox\n")
    fh.write("{moveto gsave\n")
    fh.write("dx 0 rlineto \n") 
    fh.write("0 dy rlineto \n")
    fh.write("dx neg 0 rlineto \n")
    fh.write("closepath\n")
    fh.write("0.1 setlinewitttdth\n")
    fh.write("setgray fill\n") 
    fh.write("stroke grestore}def \n")
    fh.write("\n")
    fh.write("/printMatrix\n")
    fh.write("{newpath\n")
    fh.write("0 1 rows 1 sub\n")
    fh.write("{/rowcount exch def\n")
    fh.write("0 1 columns\n")
    fh.write("{/colcount exch def\n")
    fh.write("gsave\n")
    fh.write("pagewidth colcount mul neg\n")
    fh.write("pageheight rowcount mul neg\n")
    fh.write("translate\n")

    for ix in range(numZ):
        for iy in range(numPhi):
            fh.write("/x "+ix.__str__()+" dx mul def\n")
            fh.write("/y "+iy.__str__()+" dy mul def\n")
            fh.write("/GrayScale "+GrayScale[ix][iy].__str__()+" def\n")
            fh.write("GrayScale x 40 add y TheBox\n") # x+40 for margin

    fh.write(" gsave showpage grestore grestore\n")
    fh.write("}for }for }def\n")
    fh.write("\n")
    
    fh.write("%----code-----\n")
    fh.write("\n")
    fh.write("printMatrix\n")
    fh.close()

def writepsA4(fileout,Z,R,numZ,numPhi,GrayScale):
    """
    Create and write the postscript code that print the InkShimMatrix
    in five A4 pages.
    Name of the .ps file is what had been put instead of fileout.
    """

    ### size of boxes ###
    dy=Z/float(numZ-1)             #meters
    dx=2.0*math.pi*R/float(numPhi) #meters
    dx=int(dx*72/0.0254)           #points = 1/72 inch
    dy=int(dy*72/0.0254)           #points

    ### write ps file ###
    fh=open(fileout,"w+")
    fh.write("%Postscript to print InkShimMatrix\n")
    
    fh.write("%----Variables----\n")
    fh.write("/dx "+dx.__str__()+" def\n")
    fh.write("/dy "+dy.__str__()+" def\n")
    fh.write("/inch {72 mul}def\n")
    fh.write("/pagewidth 571 def\n")        # for the margin it's not 595
    fh.write("/pageheight 842 def \n")      #A4 595*842 points
    fh.write("/rows 1 def \n")              
    fh.write("/columns 6 def \n")           #number of pages 6 
    fh.write("\n")
    
    fh.write("%----Procedures----\n")
    fh.write("\n")
    fh.write("/TheBox\n")
    fh.write("{moveto gsave\n")
    fh.write("dx 0 rlineto \n") 
    fh.write("0 dy rlineto \n")
    fh.write("dx neg 0 rlineto \n")
    fh.write("closepath\n")
    fh.write("0.1 setlinewidth\n")
    fh.write("setgray fill\n") 
    fh.write("stroke grestore}def \n")
    fh.write("\n")
    fh.write("/printMatrix\n")
    fh.write("{newpath\n")
    fh.write("0 1 rows\n")
    fh.write("{/rowcount exch def\n")
    fh.write("0 1 columns 1 sub\n")
    fh.write("{/colcount exch def\n")
    fh.write("gsave\n")
    fh.write("pagewidth colcount mul neg\n")
    fh.write("pageheight rowcount mul neg\n")
    fh.write("translate\n")

    for ix in range(numPhi): 
        for iy in range(numZ): 
            fh.write("/x "+ix.__str__()+" dx mul def\n")
            fh.write("/y "+iy.__str__()+" dy mul def\n")
            fh.write("/GrayScale "+GrayScale[iy][ix].__str__()+" def\n")
            fh.write("GrayScale x 50 add y 45 add TheBox\n") # need margins

    fh.write(" gsave showpage grestore grestore\n")
    fh.write("}for }for }def\n")
    fh.write("\n")
    
    fh.write("%----code-----\n")
    fh.write("\n")
    fh.write("printMatrix\n")
    fh.close()

if __name__== "__main__":
    writepsA3()
    writepsA4()
             

