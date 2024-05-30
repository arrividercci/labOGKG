import math
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import Button
import random
import matplotlib
matplotlib.use('TkAgg')

plt.ion()

def findDistance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

def circleBy2Points(p1, p2):
    center = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
    radius = findDistance(p1, p2) / 2
    return center, radius

def circleBy3Points(p1, p2, p3):
    A = p1[0] * (p2[1] - p3[1]) - p1[1] * (p2[0] - p3[0]) + p2[0] * p3[1] - p3[0] * p2[1]
    B = (p1[0] ** 2 + p1[1] ** 2) * (p3[1] - p2[1]) + (p2[0] ** 2 + p2[1] ** 2) * (p1[1] - p3[1]) + (p3[0] ** 2 + p3[1] ** 2) * (p2[1] - p1[1])
    C = (p1[0] ** 2 + p1[1] ** 2) * (p2[0] - p3[0]) + (p2[0] ** 2 + p2[1] ** 2) * (p3[0] - p1[0]) + (p3[0] ** 2 + p3[1] ** 2) * (p1[0] - p2[0])
    D = (p1[0] ** 2 + p1[1] ** 2) * (p3[0] * p2[1] - p2[0] * p3[1]) + (p2[0] ** 2 + p2[1] ** 2) * (p1[0] * p3[1] - p3[0] * p1[1]) + (p3[0] ** 2 + p3[1] ** 2) * (p2[0] * p1[1] - p1[0] * p2[1])
    if A == 0:
        return None, float('inf')
    center = (-B / (2 * A), -C / (2 * A))
    radius = math.sqrt(B ** 2 + C ** 2 - 4 * A * D) / (2 * abs(A))
    return center, radius

def welzl(points, boundary_points=[]):
    if len(points) == 0 or len(boundary_points) == 3:
        if len(boundary_points) == 0:
            return (0, 0), 0
        elif len(boundary_points) == 1:
            return boundary_points[0], 0
        elif len(boundary_points) == 2:
            return circleBy2Points(boundary_points[0], boundary_points[1])
        elif len(boundary_points) == 3:
            return circleBy3Points(boundary_points[0], boundary_points[1], boundary_points[2])
    p = points.pop()
    center, radius = welzl(points, boundary_points)
    if findDistance(center, p) <= radius:
        points.append(p)
        return center, radius
    boundary_points.append(p)
    center, radius = welzl(points, boundary_points)
    boundary_points.pop()
    points.append(p)
    return center, radius

fig_points, board = plt.subplots()
board.set_xlim(0, 100)
board.set_ylim(0, 100)

def drawConvexHull(points):
    xlim = board.get_xlim()
    ylim = board.get_ylim()
    if len(points) < 3:
        board.scatter(points[:, 0], points[:, 1], color='r')
        board.axis('equal')
        board.set_xlim(xlim)
        board.set_ylim(ylim)
        return

    hull = ConvexHull(points)
    for simplex in hull.simplices:
        board.plot(points[simplex, 0], points[simplex, 1], 'k-')
    board.scatter(points[:, 0], points[:, 1], color='r')
    board.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'r--', lw=2)
    board.plot([points[hull.vertices[-1], 0], points[hull.vertices[0], 0]],
               [points[hull.vertices[-1], 1], points[hull.vertices[0], 1]], 'r--', lw=2)
    board.axis('equal')
    board.set_xlim(xlim)
    board.set_ylim(ylim)

def drawCircle(center, radius):
    xlim = board.get_xlim()
    ylim = board.get_ylim()
    circle = Circle(center, radius, edgecolor='b', facecolor='none')
    board.add_patch(circle)
    board.axis('equal')
    board.set_xlim(xlim)
    board.set_ylim(ylim)

def findAndDraw(points):
    center, radius = welzl(points.tolist())
    print("Center:", center)
    print("Radius:", radius)
    drawConvexHull(points)
    drawCircle(center, radius)
    board.set_xlabel('X')
    board.set_ylabel('Y')

points = []

def onclick(event):
    x = event.xdata
    y = event.ydata
    if x is not None and y is not None:
        if event.button == 1 and event.inaxes != compute_button_ax:
            points.append((x, y))
            board.scatter(x, y, color='r')
            fig_points.canvas.draw()

def findButtonClicked(event):
    findAndDraw(np.array(points))

def randomPoints():
    xlim = board.get_xlim()
    ylim = board.get_ylim()
    board.clear()
    board.set_xlim(0, 100)
    board.set_ylim(0, 100)
    points.clear()
    randomPoints = [(random.uniform(30, 60), random.uniform(30, 60)) for _ in range(60)]
    for point in randomPoints:
        points.append(point)
        board.scatter(point[0], point[1], color='r')
    fig_points.canvas.draw()
    board.set_xlim(xlim)
    board.set_ylim(ylim)

def clearPoints():
    xlim = board.get_xlim()
    ylim = board.get_ylim()
    board.clear()
    board.set_xlim(0, 100)
    board.set_ylim(0, 100)
    points.clear()
    fig_points.canvas.draw()
    board.set_xlim(xlim)
    board.set_ylim(ylim)

clear_button_ax = fig_points.add_axes([0.2, 0.02, 0.1, 0.06])
clearButton = Button(clear_button_ax, 'Clear')
clearButton.on_clicked(lambda event: clearPoints())
random_button_ax = fig_points.add_axes([0.4, 0.02, 0.1, 0.06])
randomButton = Button(random_button_ax, 'Random')
randomButton.on_clicked(lambda event: randomPoints())
compute_button_ax = fig_points.add_axes([0.6, 0.02, 0.1, 0.06])
computeButton = Button(compute_button_ax, 'Find')
computeButton.on_clicked(findButtonClicked)
fig_points.canvas.mpl_connect('button_press_event', onclick)
plt.show()
plt.show(block=True)
