
#-------------------------------FUNCTION DEFINITION-----------------------------
def Make_a_rectangoular_crop(img_collection, topx, topy, bottomx, bottomy):

    #import matplotlib.pyplot as plt
    from PIL import Image

    img_collection_cropped = []
    #img = ((100,100))

    for i in range(len(img_collection)):
        img = img_collection[i]

        cropped  = img[topx:bottomx, topy:bottomy]
        img_collection_cropped.append(cropped)
        #plt.matshow(cropped)
        #plt.show()

    return(img_collection_cropped)


#--------------------------Contours plotting------------------------------------


def drawShape(img, coordinates, color):
        # This function draws all contours founded on the mask
        # In order to draw our line in red
    #img = skimage.color.gray2rgb(img)

        # Make sure the coordinates are expressed as integers
    #coordinates = coordinates.astype(int)

    #img[coordinates[:, 0], coordinates[:, 1]] = color

    return img


def Contours_printing(img_float, contours):
        # Display the image and plot all contours found
    from matplotlib import pyplot as plt
    import numpy as np

    fig, ax = plt.subplots()
    ax.imshow(img_float, interpolation='nearest', cmap=plt.cm.gray)

    for n, contour in enumerate(contours):
        ax.plot(contour[:, 1], contour[:, 0], linewidth=2)

    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show(block=False)

        #Create a Black mask of the same dimentions of the image
    mask = np.zeros_like(img_float) # Create mask where white is what we want, black otherwise

    #Drow contours on 'mask'
    for contour in contours:
        mask = drawShape(mask, contour, [255, 0, 0])

#------------------------Find best contour--------------------------------------

def Find_Contours(img, Contour_Limit):

    import numpy as np
    from PIL import Image
    import skimage
    from skimage import data, io, filters, measure, img_as_float, img_as_uint
    from matplotlib import pyplot as plt

    contour = []
    contour_fake = []
    max_index = 0
    #percent = 0
    #Contour_Limit = 0.20#0.15
    gap = 0.001
    flag = 0

    total = 0
    inside = 0
    # Computing the contour lines...
    contour_fake = measure.find_contours(img, Contour_Limit)
    # Select the largest contour lines
    for i in range(0,len(contour_fake)):
        if (len(contour_fake[i])> len(contour_fake[max_index])):
            max_index = i
    #print len(contour_fake)

    contour = [contour_fake[max_index]]
    appoggio_up = 0
    appoggio_down = 0
    appoggio_dx = 0
    appoggio_sx = 0
    snd = False

    for i in contour[0]:
        if i[0] == 0:
            appoggio_up = 1
        if i[0] == len(img[0])-1:
            appoggio_down = 1
        if i[1] == 0:
            appoggio_sx =1
        if i[1] == len(img[1])-1:
            appoggio_dx = 1

        if i[0]==0 or i[0]==len(img[0])-1:
            appoggio_or = 1
        if i[1]==0 or i[1]==len(img[1])-1:
            appoggio_ver = 1
    if (appoggio_down == 1 and appoggio_dx == 1) or (appoggio_up == 1 and appoggio_sx == 1) or (appoggio_down== 1 and appoggio_sx == 1) or (appoggio_up == 1 and appoggio_dx == 1):
        snd = True
    print 'snd =', snd

    #print snd

    if snd == True:
        if max_index!=0:
            snd_index = 0
        else:
            snd_index = 1

        for i in range(0,len(contour_fake)):
            if (len(contour_fake[i])>len(contour_fake[snd_index]) and i!=max_index):
                snd_index = i

        for i in range(0, len(contour_fake[snd_index])):
            contour[0] = np.append(contour[0], [contour_fake[snd_index][i]], axis = 0)

    print 'Here (Figure 2) there are contours found with Contour_Limit = %s' % str(Contour_Limit)

    return contour

#------------------------------------------------------------------------------

def Find_a_mask(img_collection, Contour_Limit):
    #This function finds and draws contours for a given image.
    #It returns the cropped image as a numpy array of float

    # Find contours the best Contour_Limit for the first image of the series
    contours = Find_Contours(img_collection[0], Contour_Limit)
    Contours_printing(img_collection[0], contours)
    return(contours)


#-------------------------------------------------------------------------------


def Inside_outside_check(point, contours):
    #Mi raccomando Contours lo si deve passare come Contours[0]

    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon

    point = Point(point)
    polygon = Polygon(contours)
    return(polygon.contains(point))


#-------------------------------------------------------------------------------

def Apply_a_mask(img_collection, Contours, DIM_X, DIM_Y):

    import numpy as np
    import matplotlib.pyplot as plt

    img_collection=np.asarray(img_collection)

    #print len(img_collection)
    #plt.matshow(img_collection)
    #plt.show()

    for x in range(DIM_X):
        for y in range(DIM_Y):
            point = (x,y)
            if not(Inside_outside_check(point, Contours[0])):
                #se il punto e' nel contorno
                img_collection[x,y] = 'nan'
    #plt.matshow(img_collection)
    #plt.show()

    return img_collection
