#!/usr/bin/env python3

'''arrange several images according to a layout'''

# import yaml
import sys
import os
import argparse
from time import sleep
from PIL import Image, ImageFont, ImageDraw, ImageChops
from collections import deque, OrderedDict
import math
font = ImageFont.load_default()


def autocrop(im):
    '''autocrop image based on color of top left pixel'''
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    else:
        return im

def addborder(im, size=5, color='white'):
    '''add border to image with given size and color'''
    if type(size) == list:
        if len(size) == 2: # left-right and top-bottom borders are given
            left, top = size
            right, bottom = size
        elif len(size) == 4: # all four borders separately given
            left, top, right, bottom = size
    else: # uniform border size given
        top = left = bottom = right = size
    if top == left == bottom == right == 0:
        return im
    w, h = im.size
    newsize = (w+left+right, h+top+bottom)
    tim = Image.new('RGBA', newsize, color)
    tim.paste(im, (left, top))
    return tim

def getconf(configfile):
    '''load configuration file into dictionary'''
    f = open(configfile)
    conf = yaml.load(f, Loader=yaml.Loader)
    print(conf)
    f.close()
    return conf

def imagelayout(conf, reportsizes=False, imagefiles=[], outputfile=''):
    '''generate layout from given configuration'''
    # initialize variables
    
    wd = conf.get('inputdir', os.getcwd())
    
    # set up image list
    imgdata = conf['images']
    
    for imgid in imgdata:
        if type(imgdata[imgid]) == str: # only filename is given
            imgdata[imgid] = {'file':imgdata[imgid]}
        if imgdata[imgid]['file'][0] == '$': # file name is taken from command line argument
            imgdata[imgid]['file'] = imagefiles[int(imgdata[imgid]['file'][1:])-1]
        if imgdata[imgid]['file'].startswith('BLANK-'): # use a blank image
            w, h = [int(s) for s in imgdata[imgid]['file'].split('-')[1].split('x')]
            im = Image.new('RGBA', (w, h), color=conf.get('paddingcolor', 'white'))
            if 'label' not in imgdata[imgid]:
                imgdata[imgid]['label'] = {'text': ''}
            elif 'text' not in imgdata[imgid]['label']:
                imgdata[imgid]['label']['text'] = ''
        else:
            im = Image.open(os.path.join(wd, imgdata[imgid]['file']))
        imgdata[imgid]['size'] = (0, 0, *im.size)
        imgdata[imgid]['parts'] = [imgid]
        imgdata[imgid]['img'] = im
    
    if reportsizes: # -s option is given, report sizes and exit
        for imgid in imgdata:
            w, h = imgdata[imgid]['size'][2:]
            print('%s %s: %dx%d aspect=%f' % (imgid, imgdata[imgid]['file'], w, h, w/h))
        sys.exit()
    
    # autocrop, crop images if requested
    
    if conf.get('autocrop', False):
        for imgid in imgdata:
            if not imgdata[imgid].get('autocrop', True):
                continue # do not autocrop this image
            imb = autocrop(imgdata[imgid]['img'])
            imgdata[imgid]['img'] = imb
            imgdata[imgid]['size'] = (0, 0, *imb.size)
    
    for imgid in imgdata:
        if imgdata[imgid].get('autocrop'):
            imb = autocrop(imgdata[imgid]['img'])
            imgdata[imgid]['img'] = imb
            imgdata[imgid]['size'] = (0, 0, *imb.size)
        if 'crop' in imgdata[imgid]: # cropping defined for individual image
            size = imgdata[imgid]['crop']
            if type(size) == list:
                if len(size) == 2:
                    left, top = size
                    right, bottom = size
                elif len(size) == 4:
                    left, top, right, bottom = size
            else:
                top = left = bottom = right = size
            if top == left == bottom == right == 0:
                continue
            w, h = imgdata[imgid]['size'][2:]
            x1, y1 = left, top
            x2, y2 = w-right, h-bottom
            im = imgdata[imgid]['img']
            imb = im.crop((x1, y1, x2, y2))
            imgdata[imgid]['img'] = imb
            imgdata[imgid]['size'] = (0, 0, *imb.size)
            
    # add borders to images if requested
    
    bordersize = 0
    bordercolor = 'white'
    if 'border' in conf:
        bordersize = conf['border'].get('size', 10)
        bordercolor = conf['border'].get('color', 'white')
    
    for imgid in imgdata:
        bs = bordersize
        bc = bordercolor
        if 'border' in imgdata[imgid]: # border defined for individual image
            bs = imgdata[imgid]['border'].get('size', bordersize)
            bc = imgdata[imgid]['border'].get('color', bordercolor)
        im = imgdata[imgid]['img']
        imb = addborder(im, bs, bc)
        imgdata[imgid]['img'] = imb
        imgdata[imgid]['size'] = (0, 0, *imb.size)
            
    # determine final sizes
    
    if 'layout' not in conf:
        if len(imgdata) == 1:
            conf['layout'] = {'hjoin': list(imgdata)}
        else:
            raise ValueError('layout not specified')
    layout = conf['layout']
    
    # topological sorting of directed graph defining the layout
    Q = deque()
    Q.appendleft(layout)
    S = [layout]
    
    c = 1
    
    while Q:
        v = Q.pop()
        if type(v) == dict:
            jointype, joinlist = list(v.items())[0]
            jointypec = jointype+str(c)
            v[jointypec] = v.pop(jointype)
            c += 1
            for w in v[jointypec]:
                if type(w) == dict and w not in S:
                    S.append(w)
                    Q.appendleft(w)
    
    adj = OrderedDict() # stores hjoin and vjoin lists
    for dic in S[::-1]:
        name = list(dic)[0]
        ilist = []
        for e in dic[name]:
            if type(e) == dict:
                ilist.append(list(e)[0])
            else:
                ilist.append(e)
        adj[name] = ilist
    
    # determine positions and sizes in final image
    
    for iid in imgdata:
        if 'fixedsize' not in imgdata[iid]:
            imgdata[iid]['fixedsize'] = False

    # do the joinings
    for name in adj:
        imgdata[name] = {'fixedsize': False}
        tomerge = adj[name]
        if True in [imgdata[iid]['fixedsize'] for iid in tomerge]:
            imgdata[name]['fixedsize'] = True
        widths = [imgdata[iid]['size'][2]-imgdata[iid]['size'][0] for iid in tomerge]
        heights = [imgdata[iid]['size'][3]-imgdata[iid]['size'][1] for iid in tomerge]
        fixs = [imgdata[iid]['fixedsize'] for iid in tomerge]
        # find largest height/width image, others will be aligned to this
        if name.startswith('vjoin'):
            im0 = widths.index(max(widths))
        elif name.startswith('hjoin'):
            im0 = heights.index(max(heights))
        # generate vertical join
        if name.startswith('vjoin'):
            wfac = [widths[im0]/widths[j] for j in range(len(widths))]
            hfac = [1 if fixs[j] else wfac[j] for j in range(len(widths))]
            newheights = [hfac[j]*heights[j] for j in range(len(widths))]
            toth = sum(newheights)
            totw = widths[im0]
            for j in range(len(tomerge)):
                kx, ky, kw, kh = imgdata[tomerge[j]]['size']
                jw = kw if fixs[j] else widths[im0]
                jh = newheights[j]
                jx = widths[im0]/2-kw/2 if fixs[j] else 0
                jy = sum(newheights[:j])
                imgdata[tomerge[j]]['size'] = (jx, jy, jw, jh)
                for part in imgdata[tomerge[j]]['parts']:
                    if part == tomerge[j]:
                        continue
                    kx, ky, kw, kh = imgdata[part]['size']
                    if fixs[j]:
                        imgdata[part]['size'] = (jx+kx, jy+ky, kw, kh)
                    else:
                        imgdata[part]['size'] = (wfac[j]*kx, jy+wfac[j]*ky, wfac[j]*kw, wfac[j]*kh)
        # generate horizontal join
        elif name.startswith('hjoin'):
            hfac = [heights[im0]/heights[j] for j in range(len(heights))]
            wfac = [1 if fixs[j] else hfac[j] for j in range(len(heights))]
            newwidths = [wfac[j]*widths[j] for j in range(len(heights))]
            toth = heights[im0]
            totw = sum(newwidths)
            for j in range(len(tomerge)):
                kx, ky, kw, kh = imgdata[tomerge[j]]['size']
                jw = newwidths[j]
                jh = kh if fixs[j] else heights[im0]
                jx = sum(newwidths[:j])
                jy = heights[im0]/2-kh/2 if fixs[j] else 0
                imgdata[tomerge[j]]['size'] = (jx, jy, jw, jh)
                for part in imgdata[tomerge[j]]['parts']:
                    if part == tomerge[j]:
                        continue
                    kx, ky, kw, kh = imgdata[part]['size']
                    if fixs[j]:
                        imgdata[part]['size'] = (jx+kx, jy+ky, kw, kh)
                    else:
                        imgdata[part]['size'] = (jx+hfac[j]*kx, hfac[j]*ky, hfac[j]*kw, hfac[j]*kh)
        # set size and parts for this join
        imgdata[name]['size'] = (0, 0, totw, toth)
        imgdata[name]['parts'] = sum((imgdata[iid]['parts'] for iid in tomerge), [])
    last = list(adj)[-1]
    totw, toth = imgdata[last]['size'][2:]
    
    # resize whole image if requested
    
    eids = [iid for iid in imgdata if not iid.startswith('hjoin') and not iid.startswith('vjoin')]
    if 'finalwidth' in conf or 'finalheight' in conf:
        if 'finalwidth' in conf:
            facw = conf['finalwidth']/totw
            if 'finalheight' in conf:
                fach = conf['finalheight']/toth
            else:
                fach = facw
        else:
            fach = conf['finalheight']/toth
            facw = fach
        for iid in eids:
            [x, y, w, h] = imgdata[iid]['size']
            imgdata[iid]['size'] = [facw*x, fach*y, facw*w, fach*h]
        totw = facw*totw
        toth = fach*toth
        
    # convert coords and sizes to integers by sum-preserving rounding
    
    for iid in eids:
        [x, y, w, h] = imgdata[iid]['size']
        imgdata[iid]['intsize'] = [round(x), round(y), round(x+w)-round(x), round(y+h)-round(y)]
    totwint = round(totw)
    tothint = round(toth)
    
    # create new image and paste each input image to its place
    
    reszmethod = conf.get('resizemethod', 'nearest')
    resmeth = {'nearest':Image.NEAREST, 'bilinear':Image.BILINEAR, 'bicubic':Image.BICUBIC,
      'lanczos':Image.LANCZOS}[reszmethod]
    
    imnew = Image.new('RGBA', (totwint, tothint), color=conf.get('paddingcolor', 'white'))
    for imgid in eids:
        [x, y, w, h] = imgdata[imgid]['intsize']
        im = imgdata[imgid]['img']
        imr = im.resize((w, h), resample=resmeth)
        imnew.paste(imr, (x, y))
    
    # add labels if requested
    
    if conf.get('labels', {'add': False}).get('add', False):
        fontname = conf['labels'].get('fontname', 'FreeSans')
        fontsize = conf['labels'].get('fontsize', 32)
        fontcolor = conf['labels'].get('fontcolor', 'black')
        defaultpos = conf['labels'].get('pos', 'top-left')
        defaultoffset = conf['labels'].get('offset', [0, 0])
        draw = ImageDraw.Draw(imnew)
        for iid in eids:
            label = imgdata[iid].get('label', iid)
            if type(label) == str:
                labeltext = label
                fn = fontname
                fs = fontsize
                fc = fontcolor
                pos = defaultpos
                offset = defaultoffset
            else:
                labeltext = imgdata[iid]['label'].get('text', iid)
                fn = imgdata[iid]['label'].get('fontname', fontname)
                fs = imgdata[iid]['label'].get('fontsize', fontsize)
                fc = imgdata[iid]['label'].get('fontcolor', fontcolor)
                pos = imgdata[iid]['label'].get('pos', defaultpos)
                offset = imgdata[iid]['label'].get('offset', defaultoffset)
            [x, y, w, h] = imgdata[iid]['intsize']
            posh, posw = pos.split('-')
            offsetx, offsety = offset
            imfont = ImageFont.truetype(font=fn, size=fs)
            if posh+posw != 'topleft':
                (textw, texth) = draw.textsize(labeltext, font=imfont)
                if posh == 'center':
                    y = y+(h-texth)/2
                elif posh == 'bottom':
                    y = y+h-texth
                if posw == 'center':
                    x = x+(w-textw)/2
                elif posw == 'right':
                    x = x+w-textw
            x += offsetx
            y += offsety
            draw.text((x, y), labeltext, fill=fc, font=imfont)
    
    # add title if requested
    
    if 'title' in conf and conf['title'].get('add', True):
        fontname = conf['title'].get('fontname', 'FreeSans')
        fontsize = conf['title'].get('fontsize', 36)
        fontcolor = conf['title'].get('fontcolor', 'black')
        bgcolor = conf['title'].get('bgcolor', 'white')
        height = conf['title'].get('height', int(fontsize*1.1))
        w, h = imnew.size
        imnew2 = Image.new('RGBA', (w, h+height), color=bgcolor)
        imnew2.paste(imnew, (0, height))
        draw = ImageDraw.Draw(imnew2)
        font = ImageFont.truetype(font=fontname, size=fontsize)
        title = conf['title']['text']
        (textw, texth) = draw.textsize(title, font=font)
        x = (w-textw)/2
        y = (height-texth)/2
        draw.text((x, y), title, fill=fontcolor, font=font)
        imnew = imnew2
        
    # add lines, if any
    
    if 'lines' in conf:
        cos30 = math.cos(15/180*math.pi)
        sin30 = math.sin(15/180*math.pi)
        width = conf['lines'].get('width', 3)
        color = conf['lines'].get('color', 'black')
        for linedic in conf['lines']['linelist']:
            coords = linedic['fromto']
            todraw = [(coords[2*i], coords[2*i+1]) for i in range(len(coords)//2)]
            lw = linedic.get('width', width)
            lc = linedic.get('color', color)
            arrowsize = linedic.get('arrowsize', 0)
            draw = ImageDraw.Draw(imnew)
            if arrowsize > 0:
                # draw an arrowhead
                x1, y1, x2, y2 = coords[-4:]
                vx, vy = x1-x2, y1-y2
                v = math.sqrt(vx**2+vy**2)
                vx, vy = arrowsize*vx/v, arrowsize*vy/v
                vxL = cos30*vx-sin30*vy
                vyL = sin30*vx+cos30*vy
                vxR = cos30*vx+sin30*vy
                vyR = -sin30*vx+cos30*vy
                todraw += [(x2+vxL, y2+vyL), (x2, y2), (x2+vxR, y2+vyR), (x2, y2)]
            draw.line(todraw, fill=lc, width=lw, joint='curve')
                
    # add global labels, if any
    
    if 'globallabels' in conf:
        fontname = conf['globallabels'].get('fontname', 'FreeSans')
        fontsize = conf['globallabels'].get('fontsize', 36)
        fontcolor = conf['globallabels'].get('fontcolor', 'black')
        for labdic in conf['globallabels']['labellist']:
            fn = labdic.get('fontname', fontname)
            fs = labdic.get('fontsize', fontsize)
            fc = labdic.get('fontcolor', fontcolor)
            al = labdic.get('align', 'left')
            [x, y] = labdic.get('coords', [0, 0])
            draw = ImageDraw.Draw(imnew)
            font = ImageFont.truetype(font=fn, size=fs)
            draw.text((x, y), labdic['text'], fill=fc, font=font, align=al)
    
    # add global border if requested
    
    if 'globalborder' in conf:
        imnew = addborder(imnew, size=conf['globalborder'].get('size', 20),
          color=conf['globalborder'].get('color', 'white'))
    
    # save result
    if conf['outputfile'][0] == '$':
        conf['outputfile'] = outputfile # from command line argument
    ofile = conf['outputfile'].lower()
    background = Image.new('RGBA', imnew.size, color=conf.get('paddingcolor', 'white'))
    imnew = Image.alpha_composite(background, imnew)
    # pixelscaling as last step before saving
    if 'pixelscaling' in conf:
        if type(conf['pixelscaling']) == list:
            xscale, yscale = conf['pixelscaling']
        else:
            xscale = yscale = conf['pixelscaling']
        h, w = imnew.size
        imnew = imnew.resize((int(xscale*h), int(yscale*w)), resample=resmeth)
    if ofile.endswith('.tif') or ofile.endswith('.tiff'):
        imnew.save(conf['outputfile'], compression='tiff_adobe_deflate')
    else:
        imnew.save(conf['outputfile'])
    print('Output written to %s' % (conf['outputfile']))


# if __name__ == '__main__':
#     config = {
#     "inputdir": "/Users/sanjay/Monash/Master_thesis/lab_work/Li_Lab/cluster/HLA-PepClust/data/ref_data/Gibbs_motifs_human/",
#     "images": {
#         "A": {
#             "file": "HLA_A0101.png",
#             "label": {
#                 "text": "A"
#             }
#         },
#         "B": {
#             "file": "HLA_A0201.png",
#             "label": {
#                 "text": "H",
#                 "fontsize": 50,
#                 "fontcolor": "red",
#                 "pos": "center-left"
#             }
#         },
#         "C": "HLA_A0203.png",
#         "D": "HLA_A0301.png",
#         "E": "HLA_A1102.png",
#         "F": "HLA_A0211.png",
#         "G": "BLANK-640x480"
#     },
#     "layout": {
#         "vjoin": [
#             {"hjoin": ["A", "B", "C"]},
#             {"hjoin": ["G", "G", "G"]}
#         ]
#     },
#     "resizemethod": "nearest",
#     "outputfile": "output_imggrid_dict2.png",
#     "finalwidth": 1024,
#     "finalheight": 800
# }
#     imagelayout(config, reportsizes=False, imagefiles=[], outputfile='output_imggrid_dict2.png')