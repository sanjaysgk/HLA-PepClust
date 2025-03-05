
## Main program

def main():
    
    config = {
    "inputdir": "/Users/sanjay/Monash/Master_thesis/lab_work/Li_Lab/cluster/HLA-PepClust/data/ref_data/Gibbs_motifs_human/",
    "images": {
        "A": {
            "file": "HLA_A0101.png",
            "label": {
                "text": "A"
            }
        },
        "B": {
            "file": "HLA_A0201.png",
            "label": {
                "text": "H",
                "fontsize": 50,
                "fontcolor": "red",
                "pos": "center-left"
            }
        },
        "C": "HLA_A0203.png",
        "D": "HLA_A0301.png",
        "E": "HLA_A1102.png",
        "F": "HLA_A0211.png",
        "G": "BLANK-640x480"
    },
    "layout": {
        "vjoin": [
            {"hjoin": ["A", "B", "C"]},
            {"hjoin": ["G", "G", "G"]}
        ]
    },
    "resizemethod": "nearest",
    "outputfile": "output_imggrid_dict.png",
    "finalwidth": 1024,
    "finalheight": 800
}
    # parser = argparse.ArgumentParser(description='Arrange images into a layout',
    #   epilog='See man page at https://github.com/aszilagyi/imagelayout/docs/manu.md')
    # parser.add_argument('-w', '--watch', action='store_true',
    #   help='Watch the config file and re-run upon detecting a change')
    # parser.add_argument('-s', action='store_true', help='report image sizes and exit')
    # parser.add_argument('-o', dest='outputfile', help='Output image file name (optional)')
    # parser.add_argument('configfile', help='config file')
    # parser.add_argument('imagefile', nargs='*', help='image files (optional)')
    # a = parser.parse_args()

    # imagelayout(getconf(a.configfile), reportsizes=a.s, imagefiles=a.imagefile, 
    #   outputfile=a.outputfile)
    imagelayout(config, reportsizes=False, imagefiles=[], outputfile='output_imggrid_dict.png')
    
    
    # if a.watch:
    #     print('Watching config file (%s) for changes, press Ctrl-C to quit...' % (a.configfile))
    #     mtime = os.stat(a.configfile).st_mtime
    #     while True:
    #         sleep(0.5)
    #         mt = os.stat(a.configfile).st_mtime
    #         if mt != mtime:
    #             print('Re-running...')
    #             try:
    #                 imagelayout(getconf(a.configfile), reportsizes=a.s, imagefiles=a.imagefile,
    #                   outputfile=a.outputfile)
    #             except ValueError as e:
    #                 print('ValueError:', e)
    #             mtime = mt

# if __name__ == '__main__':
#     main()