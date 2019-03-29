import os
import xml.etree.ElementTree as ET

def unzip():
    os.system('unzip psas_rocket.ork')

def zipit():
    os.system('zip -m psas_rocket.ork rocket.ork')

def body_r(eng_r):
    b_d = eng_r*2 + 0.014
    return b_d/2 #* 1000 #mm

def update_body(eng_r, l_o, l_f):
    #convert lengths to mm
    #l_o *= 1000
    #l_f *= 1000
    unzip()
    tree = ET.parse('rocket.ork')
    root = tree.getroot()
    for child in root.iter():
        #set radius
        for kid in child.iterfind('aftradius'):
            kid.text = str(body_r(eng_r))
        for kid in child.iterfind("*[name='LOX Tank']"):
            kid.find('length').text = str(l_o)
        for kid in child.iterfind("*[name='Fuel Tank']"):
            kid.find('length').text = str(l_f)
    tree.write('rocket.ork')
    zipit()
    print("Rocket modified!")
    
