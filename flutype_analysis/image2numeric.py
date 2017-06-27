"""
Module for spot quantification from images.

Also implement quality control of given spot and array/image.

"""
import matplotlib.pyplot as plt
import cv2
import numpy as np
from scipy import spatial
from analysis import base, utils
from scipy.stats import norm


def draw_imag(imag, name, **kwargs):
    fig, ax = plt.subplots(figsize=(30, 10))
    plt.title(name)
    plt.imshow(imag, **kwargs)
    plt.colorbar()
    plt.show()

def get_grid(shape, pitch, center_x, center_y, rotation=0):
    """
    returns x,y coordinates of a grid with shape, pitch, center x , center y, rotation
    """
    x_spots, y_spots = np.meshgrid(
        (np.arange(shape[1]) - (shape[1] - 1) / 2.) * pitch,
        (np.arange(shape[0]) - (shape[0] - 1) / 2.) * pitch, indexing='ij')
    theta = rotation / 180. * np.pi
    x_spots = x_spots * np.cos(theta) - y_spots * np.sin(theta) + center_x
    y_spots = x_spots * np.sin(theta) + y_spots * np.cos(theta) + center_y
    return x_spots, y_spots

def get_nearest_on_grid(self,x_grid, y_grid, xe, ye, max=float('inf')):
    """
    returns x and y coordinates and index of the nearest spot on grid to a points with (xe,ye) coordinates. xe,ye can be an array of points.
    """

    xy_grid = np.vstack([x_grid, y_grid]).reshape(2, -1).T
    xy_circles = np.vstack([xe, ye]).reshape(2, -1).T
    index_nearest = self.nearest_neighbour(xy_circles, xy_grid)

    # dist = np.linalg.norm(np.array([xy_grid[index_nearest,0]-xe,xy_grid[index_nearest,1]-ye]))

    return xy_grid[index_nearest, 0], xy_grid[index_nearest, 1], index_nearest


def nearest_neighbour(points_a, points_b):
    """
    returns coordinates of closest point in array (points a) to every point in array (points b)
    """
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[1]


def get_mean_distance(x1, y1, x2, y2):
    """
    mean dinstance bestween (x1,y1) and (x2,y2)
    """
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2).mean()


def err_func(self,params, xe, ye, gridshape):
    """
    params: Grid parameters

    returns the mean distance between leuchtefix spots in images and on grid.

    """
    pitch, center_x, center_y, rotation = params
    x_grid, y_grid = self.get_grid(gridshape, pitch, center_x, center_y, rotation)
    x_spots, y_spots, _ = self.get_nearest_on_grid(x_grid, y_grid, xe, ye)
    return self.get_mean_distance(x_spots, y_spots, xe, ye)


def err_func2(params, pitch, xe, ye, index_grid, gridshape):
    """
    params: Grid parameters

    returns the mean distance between leuchtefix spots in images and on grid.

    """
    center_x, center_y, rotation = params
    x_grid, y_grid = get_grid(gridshape, pitch, center_x, center_y, rotation)
    pts = np.vstack([x_grid, y_grid]).reshape(2, -1).T
    return get_mean_distance(pts[index_grid, 0], pts[index_grid, 1], xe, ye)



def detect_circles(imag):
    """
    detects circles from images
    """

    gray_blur = cv2.medianBlur(imag, 29)  # Remove noise before laplacian

    gray_contr, contours = get_contour(gray_blur)

    gray_lap = cv2.Laplacian(gray_contr, cv2.CV_8UC1, ksize=5)

    dilate_lap = cv2.dilate(gray_lap, (3, 3))

    # Fill in gaps from blurring. This helps to detect circles with broken edges.
    # Furture remove noise introduced by laplacian. This removes false pos in space between the two groups of circles.


    lap_blur = cv2.bilateralFilter(dilate_lap, 5, 9, 9)

    circles = cv2.HoughCircles(lap_blur, cv2.HOUGH_GRADIENT, 1, 20,
                               param1=50, param2=30, minRadius=0,
                               maxRadius=100)  # cimg = draw_circles(gray, circles)    #cimg = draw_circles(image_path, circles)
    print("{} circles detected.".format(circles[0].shape[0]))
    return circles,lap_blur


def get_contour(imag):
    """
    finds contours in image
    """
    ret, thresh = cv2.threshold( imag, 60, 255, 1)
    im2, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    return im2,contours

def draw_imag(imag ,name, **kwargs):
    fig, ax = plt.subplots(figsize=(30,10))
    plt.title(name)
    plt.imshow(imag, **kwargs)
    plt.colorbar()
    plt.show()

def draw_circles(imag, circles, **kwargs):
    """
    draws a circles in an image (img)
    """
    plt.subplots(figsize=(30, 10))
    plt.imshow(imag, **kwargs)
    plt.colorbar()
    plt.scatter(circles[0, :, 1], circles[0, :, 0], s=circles[0, :, 2], color="r")
    plt.show()

def draw_grid(imag, X, Y, name, **kwargs):

    """
    draws a circles in an image (img)
    """

    fig, ax = plt.subplots(figsize=(30, 10))
    plt.imshow(imag, **kwargs)
    plt.colorbar()
    plt.title(name)
    plt.scatter(X, Y, alpha=0.5, color="r")
    plt.show()


def spot_quality_shape(x,y):
    x_spread=(x.max()-x.min())
    y_spread=(y.max()-y.min())

    #X_hist, bin_edges_Y = np.histogram(X, density=True)

    #(mu, sigma) = norm.fit(X_hist)


    #Y_hist, bin_edges_Y = np.histogram(Y, density=True)




    return x_spread/y_spread


def spots_close_to_grid(x, y, pitch):
    near_pts = []

    ptss = np.vstack([x, y]).reshape(2, -1).T
    for i in range(len(x)):
        pts = ptss
        point = pts[i]
        pts = np.delete(pts, i, 0)

        near =nearest_neighbour(point, pts)
        if get_mean_distance(pts[near, 0], pts[near, 1], point[0], point[1]) < 3 * pitch:
            near_pts.append(i)
    pts = ptss[near_pts]

    return pts[:, 0], pts[:, 1]


class Image2numeric(base.Base):

    def  __init__(self , d_data):
        base.Base.__init__(self , d_data)
        self.imag = d_data["tif"]










    '''
    @staticmethod
    def inverte(imagem):
        imagem = (255 - imagem)
        return (imagem)
    
    '''
    @staticmethod
    def inverte(imagem):
        imagem = (255 - imagem)
        return (imagem)







