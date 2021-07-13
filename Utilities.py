import numpy as np

def SaveArrayTXT(filename='DefaultArrayName.txt', array=None):
    if array is not None:
        np.savetxt(filename,array,fmt='%d',delimiter='\t')

def ArrayToXML(array):
    import xml.etree.ElementTree as ET
    from xml.dom import minidom
    nameProtocol = input('Input the name of the protocol :').strip()
    descriptionProtocol = input('Input a short description of the protocol:').strip()
    spreadFileName = input('Input the corresponding spread file :').strip()
    protocol = ET.Element('Protocol')
    name = ET.SubElement(protocol,'Name')
    name.text = nameProtocol
    description = ET.SubElement(protocol,'Description')
    description.text = descriptionProtocol
    arraycode = ET.SubElement(protocol,'Arraycode')
    arraycode.text = '15'
    spreadfile = ET.SubElement(protocol,'SpreadFile')
    spreadfile.text = spreadFileName
    sequence = ET.SubElement(protocol,'Sequence')
    ABInj = np.unique(array[:,:2],axis=0)
    for AB in ABInj:
        measure = ET.SubElement(sequence,'Measure')
        A, B = AB[0], AB[1]
        tx = ET.SubElement(measure, 'Tx')
        tx.text = f'{A} {B}'
        idxCurr = np.where((array[:,:2] == (A, B)).all(axis=1))
        idxCurr = np.asarray(idxCurr)
        for i in idxCurr.flatten():
            rx = ET.SubElement(measure,'Rx')
            M = array[i,2]
            N = array[i,3]
            rx.text = f'{M} {N}'
    roughXML = ET.tostring(protocol)
    niceXML = minidom.parseString(roughXML)
    niceXML = niceXML.toprettyxml(indent="\t")
    fileName = nameProtocol.title().replace(' ', '') + '.xml'
    file = open(fileName, "w")
    file.write(niceXML)
    
def SortArray(array=None, type='ABMN'):
    if array is not None:
        if type == 'ABMN':
            array = array[array[:,3].argsort()]
            array = array[array[:,2].argsort(kind='mergesort')]
            array = array[array[:,1].argsort(kind='mergesort')]
            array = array[array[:,0].argsort(kind='mergesort')]
        if type == 'NMBA':
            array = array[array[:,0].argsort()]
            array = array[array[:,1].argsort(kind='mergesort')]
            array = array[array[:,2].argsort(kind='mergesort')]
            array = array[array[:,3].argsort(kind='mergesort')]
    return array

def CreateReciprocals(array):
    '''Create the reciprocals array for a given array in input.
    '''
    if array is not None:
        return SortArray(array[:,[2, 3, 0, 1]])

def SampleReciprocals(array, sampling=0.1):
    '''Sample the reciprocals array for the injections
    '''
    if array is not None:
        reciprocals = CreateReciprocals(array)
        ABInj = np.unique(reciprocals[:,:2], axis=0)
        nbSamples = int(len(ABInj)*sampling)
        print(f'{nbSamples} injection in reciprocals')
        idxKeep = np.random.choice(np.arange(len(ABInj)),nbSamples,replace=False)
        ABInjKeep = ABInj[idxKeep,:]
        ReciArray = np.asarray([0, 0, 0, 0])
        for AB in ABInjKeep:
            # Find the corresponding arrays:
            A, B = AB[0], AB[1]
            idxCurr = np.where((reciprocals[:,:2] == (A, B)).all(axis=1))
            ReciArray = np.vstack([ReciArray, np.squeeze(reciprocals[idxCurr,:])])
        ReciArray = ReciArray[1:,:]
        return SortArray(ReciArray)
    