import cmath
import math
import random

# change window size to your needs!
screen_x = 1280
screen_y = 720


class P_disc:
    """
    A class of static methods for calculations on the pioncaré disc
    """

    @staticmethod
    def translation(z, a):
        """
        input: 
            z = (complex) Point to be translated
            a = (complex) translation vector
        output:  translated Point(complex) z'
        """
        return (z + a) / (1 + a.conjugate() * z)

    @staticmethod
    def intuitive_translation_constructor(a, b):
        """
        The goal of this function is to provide an intuitive translation,
        where the first selected point (a) will be translated onto the second point (b) that was selected.
        :param a: the first translated point
        :param b: second selected point
        :return: translation function
        """

        def intuitive_translation(z):
            """
            this is a constructed method, that translates a certain given point by a predefined vector
            :param z: the point to be translated
            :return: complex translated point
            """
            first_trans = (z - a) / (1 - a.conjugate() * z)  # translate a -> 0
            return (first_trans + b) / (1 + b.conjugate() * first_trans)  # translate 0 -> b

        return intuitive_translation

    @staticmethod
    def dist(a, b):
        """
        input :
            a,b = Point(complex)
        output: (float)  hyperbolic distance between these two points on the poincaré disc
        """
        delta = abs((a - b) / (1 - b.conjugate() * a))
        return 1 / 2 * log((1 + delta) / (1 - delta))

    @staticmethod
    def equilateral_Shape(n, r, circles=False):
        """
        input :
            n = (integer) amount of Edges of the equilateral shape
            r = (float) radius of the cicles to calculate the shape corners
            circles = (boolean) debugging parameter to draw the cicles which are used to calcuated shape corners
        output: (instance of Hyperbolic_Shape) n-shape with equilateral properties on the poincare disc
        """
        global radius
        alpha = radians(float(360 / n))
        if circles == True:  # debug use , prints cicles to calcuated shape corners
            i = 1
            stroke(0, 0, 0)
            noFill()
            while i <= n:
                p = interpret_point(complex(cos(i * alpha), sin(i * alpha)) * complex(0, 1))
                i += 1
                point(p[0], p[1])
                circle(p[0], p[1], r * 2 * radius)

        d = (abs(complex(cos(2 * alpha), sin(2 * alpha))
                 - complex(cos(alpha), sin(alpha))))
        if d <= 2 * r:  # early checks if cicles intersect
            intersections = []
            i = 1
            while i <= n:  # calculated intersection with cicle_intersections() for each circlepair
                solutions = cicle_intersections(complex(cos((i + 1) * alpha), sin((i + 1) * alpha)), r,
                                                complex(cos(i * alpha), sin(i * alpha)), r, d)
                if abs(solutions[0]) <= abs(solutions[1]) <= 1:  # check which solution is closer to the center
                    intersections.append(
                        solutions[0] * complex(0, 1))  # add rotation of 90 degree (gives me inner peace this way)
                elif abs(solutions[1]) <= 1:
                    intersections.append(solutions[1] * complex(0, 1))
                else:  # None of the solutions are on the unit cicle. Will be true for all pairs, so interupt the whole loop
                    break
                i += 1
            if intersections == []:
                textbox_print("Radius to big!")
                return None
            k = len(intersections) - 1
            equilateral = (Hyperbolic_Shape(
                Hyperbolic_Line(intersections[0], intersections[k])))
            colorlist = [[0, 0, 255], [0, 255, 0], [255, 0, 0], [204, 204, 0], [157, 0, 157]]
            while k >= 1:
                picker = k % len(colorlist)
                (equilateral.include(Hyperbolic_Line(
                    intersections[k], intersections[k - 1], colorlist[picker])))
                k -= 1
            redraw_all()
            return equilateral
        else:
            textbox_print("Radius to small!")
            return None


class HyperbolicPlaneWidget:
    """
    Class of the Plane widget, combining the plane and its controls into one easily movable object
    """
    def __init__(self, x_coor=50, y_coor=50, size=(500, 280), widget_color=(100, 0, 200)):
        self.position = (x_coor, y_coor)
        self.size = size
        self.button_size = min(self.size[0]//12, self.size[1]//5)
        self.widget_identifying_color = widget_color

        x_coor += self.size[0]//12
        self.plane = HyperbolicPlane(x_coor, self.position[1] + self.button_size,
                                     size=(4*self.size[0]//5, 4*self.size[1]//5),
                                     widget_ident_color=self.widget_identifying_color)

        # init buttons
        self.button_drag_drop = DragButton(x_coor - self.button_size, y_coor, self.button_size,
                                           img_blank, img_blank_pressed)
        self.button_zoom_in = UniversalSimpleButton(x_coor, y_coor, self.button_size, plus_image, plus_pressed_image,
                                                    self.plane.change_scope_zoom_multiplicative, False, (2, 2))
        self.button_zoom_out = UniversalSimpleButton(x_coor + self.button_size, y_coor, self.button_size,
                                                     minus_image, minus_pressed_image,
                                                     self.plane.change_scope_zoom_multiplicative, False, (0.5, 0.5))
        self.button_move_up = UniversalSimpleButton(x_coor + self.button_size*2, y_coor, self.button_size,
                                                    up_image, up_pressed_image,
                                                    self.plane.move_scope, False, (0, -1))
        self.button_move_down = UniversalSimpleButton(x_coor + self.button_size*3, y_coor, self.button_size,
                                                      down_image, down_pressed_image,
                                                      self.plane.move_scope, False, (0, 1))
        self.button_move_right = UniversalSimpleButton(x_coor + self.button_size*4, y_coor, self.button_size,
                                                       right_image, right_pressed_image,
                                                       self.plane.move_scope, False, (1, 0))
        self.button_move_left = UniversalSimpleButton(x_coor + self.button_size*5, y_coor, self.button_size,
                                                      left_image, left_pressed_image,
                                                      self.plane.move_scope, False, (-1, 0))
        self.button_rotate_counter_clock = UniversalSimpleButton(x_coor + self.button_size*6, y_coor, self.button_size,
                                                      rotate_counter_clock, rotate_counter_clock,
                                                      self.plane.change_orientation, False, QUARTER_PI/2)

        self.button_axis_switch = UniversalSimpleButton(x_coor + self.button_size * 7, y_coor,
                                                        self.button_size,
                                                        img_axis, img_axis_pressed,
                                                        self.switch_axis, True, 0)
        self.button_coor_switch = UniversalSimpleButton(x_coor + self.button_size * 8, y_coor,
                                                        self.button_size,
                                                        img_coordinate, img_coordinate_pressed,
                                                        self.switch_axis_steps, True, 0)
        self.button_grid_switch = UniversalSimpleButton(x_coor + self.button_size * 9, y_coor,
                                                        self.button_size,
                                                        img_grid, img_grid_pressed,
                                                        self.switch_grid, True, 0)

        self.button_list = [self.button_move_down, self.button_move_up, self.button_move_right, self.button_move_left,
                            self.button_zoom_in, self.button_zoom_out, self.button_rotate_counter_clock,
                            self.button_axis_switch, self.button_coor_switch, self.button_grid_switch,
                            self.button_drag_drop]

        plane_list.append(self)
        self.draw()

    def draw_border(self):
        """
        draws the border of the widget
        :return: None
        """
        stroke(self.widget_identifying_color[0], self.widget_identifying_color[1], self.widget_identifying_color[2])
        fill(0, 0, 0)
        rectMode(CORNER)
        rect(self.position[0], self.position[1], self.size[0], self.size[1])

    def draw(self):
        """
        draws the initial plane of the widget
        :return: None
        """
        self.draw_border()
        for button in self.button_list:
            button.update()

        self.plane.draw_plane()

    def redraw(self, selector=-1):
        """
        redraws the plane of the widget
        :param selector: (int) if given redraws only 1 state, else redraws all states
        :return: None
        """
        self.plane.redraw_plane(selector)

    def redraw_all(self, selector=-1):
        """
        redraws the whole widget
        :param selector: (int) if given redraws only 1 state, else redraws all states
        :return:
        """
        self.draw_border()
        for button in self.button_list:
            button.update()
        self.redraw(selector)

    def switch_axis(self):
        self.plane.show_axis = not self.plane.show_axis
        self.redraw_all()

    def switch_axis_steps(self):
        self.plane.show_axis_steps = not self.plane.show_axis_steps
        self.redraw_all()

    def switch_grid(self):
        self.plane.show_grid = not self.plane.show_grid
        self.redraw_all()

    def drag(self, vector):
        """
        changes the position of every element in the widget an of the widget itself by a vector
        :param vector: (int, int)
        :return: None
        """
        self.position = (self.position[0] + vector[0], self.position[1] + vector[1])
        self.plane.plane_position = (self.plane.plane_position[0] + vector[0], self.plane.plane_position[1] + vector[1])
        for button in self.button_list:
            button.px += vector[0]
            button.py += vector[1]
        redraw_all()


class HyperbolicPlane:
    """
    Class of the hyperbolic plane scope, capable of traversing the hyperbolic plane.
    """

    def __init__(self, x_coor= screen_x // 2, y_coor=screen_y // 2, size=(400, 200), widget_ident_color=(100, 0, 200)):
        # properties:
        self.plane_position = (x_coor, y_coor)  # drawing position of the plane
        self._scope_position = (0.0, 0.0)  # the position of the scope to move through the hyperbolic plane
        self.scope = size  # the size of the view, since the whole hyperbolic plane cant be displayed, todo: rename to plane_size
        self._scope_zoom = (1.0, 1.0)  # the zooming factor of the scope, be careful with useage, since it might falsify data (rounding errors)
        self._orientation = 0.0  # orientation of the scope towards the circle

        self._x_step_dist = 20  # step distance for the euclidian grid on the x-axis (int)
        self._y_step_dist = 20  # step distance for the euclidian grid on the y-axis (int)
        self._set_step_distance()

        # borders of scope
        self._min_x_value = -10.0
        self._min_y_value = 0.0
        self._max_x_value = 10.0
        self._max_y_value = 10.0
        self._set_border()

        # settings:
        self.show_grid = False  # show an euclidian grid, not supported atm
        self.show_axis = False  # show the x- and y-axis
        self.show_axis_steps = False  # show euclidian marks on the x- and y-axis

        self.background_color = (255, 255, 255)
        self.axis_color = (0, 0, 0)
        self.grid_color = (0, 0, 0)
        self._ident_color = widget_ident_color  # color to identify the plane by

    def draw_plane(self):
        """
        drawing the hyperbolic plane as a rectangle
        :return: None
        """
        noStroke()
        fill(self.background_color[0], self.background_color[1], self.background_color[2])
        rectMode(CORNER)
        rect(self.plane_position[0], self.plane_position[1], self.scope[0], self.scope[1])

        if self.show_grid:
            self.draw_grid()
        if self.show_axis:
            self.draw_coordinates()

        self._draw_orientation()

    def redraw_plane(self, selector=-1):
        """
        function to redraw the plane, the difference to the normal draw function is,
        that the hyperbolic lines also get redrawn
        :selector: if a selector is given, only lines in that selector will be redrawn TODO:
        :return: None
        """
        self.draw_plane()

        if selector >= 0:
            if len(shapelist) > selector:  # check if selected shape exist for error catching
                for stored_line in shapelist[selector].lines:
                    stored_line.drawline_on_plane(self)
            return None

        # draw lines on plane
        for shape in shapelist:
            for stored_line in shape.lines:
                stored_line.drawline_on_plane(self)

    def draw_coordinates(self):
        """
        drawing axis with euclidian markings onto the hyperbolic plane
        :return: None
        """
        noFill()
        stroke(self.axis_color[0], self.axis_color[1], self.axis_color[2])

        # calc origin point
        x_zero_coor = self.plane_position[0] + self.scope[0] // 2 - self._scope_position[0] * self._x_step_dist
        y_zero_coor = self.plane_position[1] + self.scope[1] - self._scope_position[1] * self._y_step_dist

        # calc bottom left first (int, int) absolute position
        down_left_start_point = self.get_point_relative_to_plane(complex(ceil(self._min_x_value),
                                                                         ceil(self._min_y_value)))

        # drawing the x-axis:
        #  note: since the x-axis is always at the bottom, and you cant move past it,
        #  this first case is essentially the same as the second one,
        #  but just in case someone wants to fiddle with the movement I will leave the distinction in here
        if self._min_y_value <= 0 <= self._max_y_value:  # checking if x-axis is in scope
            line(self.plane_position[0], y_zero_coor,
                 self.plane_position[0] + self.scope[0], y_zero_coor)

        else:  # x-axis moved out of scope, draw visual cues at bottom
            y_zero_coor = self.plane_position[1] + self.scope[1]

        # drawing the y-axis:
        if self._min_x_value <= 0 <= self._max_x_value:  # checking if y-axis is in scope
            line(x_zero_coor, self.plane_position[1] + self.scope[1],
                 x_zero_coor, self.plane_position[1])

        else:  # y-axis out of scope, draw visual cues on the side, the y-axis is at
            x_zero_coor = self.plane_position[0] + (self._min_x_value < 0)*self.scope[0]

        if self.show_axis_steps:

            # drawing visual cues for x-axis
            for step in range(0, self.scope[0], self._x_step_dist):
                line(down_left_start_point[0] + step,
                     y_zero_coor,
                     down_left_start_point[0] + step,
                     y_zero_coor - self.scope[1] / 25)

            # drawing visual cues for y-axis
            for step in range(0, self.scope[1], self._y_step_dist):
                line(x_zero_coor - self.scope[0] / 50,
                     down_left_start_point[1] - step,
                     x_zero_coor + self.scope[0] / 50,
                     down_left_start_point[1] - step)

    def draw_grid(self):
        """
        drawing an euclidian grid onto the hyperbolic plane, not working atm
        :return: None
        """
        noFill()
        stroke(self.grid_color[0], self.grid_color[1], self.grid_color[2])
        # x-axis
        for step in range(0, self.scope[0], self._x_step_dist):
            line(self.plane_position[0] + step,
                 self.plane_position[1] + self.scope[1],
                 self.plane_position[0] + step,
                 self.plane_position[1])

        # y-axis
        for step in range(0, self.scope[1], self._y_step_dist):
            line(self.plane_position[0],
                 self.plane_position[1] + self.scope[1] - step,
                 self.plane_position[0] + self.scope[0],
                 self.plane_position[1] + self.scope[1] - step)

    def switch_axis(self):
        """
        Switches the truth value of show_axis
        :return: None
        """
        self.show_axis = not self.show_axis

    def switch_axis_steps(self):
        """
        Switches the truth value of show_axis_steps
        :return: None
        """
        self.show_axis_steps = not self.show_axis_steps

    def switch_grid(self):
        """
        Switches the truth value of show_grid
        :return: None
        """
        self.show_grid = not self.show_grid

    def _draw_orientation(self):
        """
        drawing the orientation onto the circle by coloring a line at the outside of the circle,
        the lenght of the marker is hardcoded
        :return: None
        """
        stroke(self._ident_color[0], self._ident_color[1], self._ident_color[2])
        noFill()
        d = 20  # length of marker
        a = (origin[0] + math.cos(self._orientation)*(radius+5), origin[1] - math.sin(self._orientation)*(radius+5))
        b = (origin[0] + math.cos(self._orientation)*(radius+d), origin[1] - math.sin(self._orientation)*(radius+d))

        line(a[0], a[1], b[0], b[1])

    def get_point_relative_to_plane(self, point):
        """
        get the absolute position of a point in the app space, based on the relative coordinates in the plane
        :param point: relative point in plane: complex
        :return: absolute point in app tuple[int, int]
        """
        absolute_x_pos = round(self.plane_position[0]
                               - self._scope_position[0]*self._x_step_dist
                               + self.scope[0]/2 + point.real * self._x_step_dist)
        absolute_y_pos = round(self.plane_position[1] - self._scope_position[1]*self._y_step_dist +
                               self.scope[1] - point.imag * self._y_step_dist)
        return absolute_x_pos, absolute_y_pos

    def get_scope_relative_scaling(self, plane_vector):
        """
        scales a vector (e.g. distance) to fit the scaling of the plane
        :param plane_vector: the vector that needs to be adjusted
        :return: a vector (tuple[float, float]) scaled in relation to the scope
        """
        scaled_x_coor = plane_vector.real * self._x_step_dist
        scaled_y_coor = plane_vector.imag * self._y_step_dist
        return scaled_x_coor, scaled_y_coor

    def _set_border(self):
        """
        calculates the border values of the scope
        :return: None
        """
        # borders of scope
        self._min_x_value = self._scope_position[0] - float(self.scope[0]) / float(self._x_step_dist*2)
        self._min_y_value = -self._scope_position[1]
        self._max_x_value = self._scope_position[0] + float(self.scope[0]) / float(self._x_step_dist*2)
        self._max_y_value = -self._scope_position[1] + float(self.scope[1]) / float(self._y_step_dist)

    def _set_step_distance(self):
        """
        calculates the step distance in pixel, calcs as zoom times 20, meaning default is 1 step 20 pixel
        :return: None
        """
        self._x_step_dist = int(20 * self._scope_zoom[0])  # step distance for the euclidian grid on the x-axis
        self._y_step_dist = int(20 * self._scope_zoom[1])  # step distance for the euclidian grid on the y-axis

    def set_orientation(self, angle):
        """
        setter for the orientation
        :param angle: float value from -PI to PI
        :return: None
        """
        if 0 <= angle < TWO_PI:
            self._orientation = angle
        redraw_all(selector)

    def change_orientation(self, angle):
        """
        changes the orientation additive
        :param angle: amount by which to change the orientation radians
        :return: None
        """
        self.set_orientation((self._orientation + angle) % TWO_PI)

    def get_orientation(self):
        """
        getter of plane orientation
        :return:
        """
        return self._orientation

    def get_scope_zoom(self):
        return self._scope_zoom

    def change_scope_zoom_additive(self, zoom):
        """
        changes the zoom of the plane by the amount given, resulting zoom cant be smaller than 0
        :param zoom: tuple[float, float] to set the zoom of the plane
        :return: None
        """
        self.set_scope_zoom((self._scope_zoom[0] + zoom[0], self._scope_zoom[1] + zoom[1]))

    def change_scope_zoom_multiplicative(self, zoom_multiplier):
        """
        changes the zoom of the scope by a factor
        :param zoom_multiplier: tuple[float, float] both float > 0
        :return: None
        """
        self.set_scope_zoom((self._scope_zoom[0] * zoom_multiplier[0], self._scope_zoom[1] * zoom_multiplier[1]))

    def set_scope_zoom(self, zoom):
        """
        changes the zoom of the scope to a different zoom
        :param zoom: the new zoom tuple[float, float], both float > 0
        :return: None
        """
        if 2^54 < zoom[0] < 0.05:
            zoom = (self._scope_zoom[0], zoom[1])
        if 2^54 < zoom[1] < 0.05:
            zoom = (zoom[0], self._scope_zoom[1])
        self._scope_zoom = zoom

        # x,y step distance
        self._set_step_distance()

        # borders of scope
        self._set_border()

        self.redraw_plane(selector)

    def move_scope(self, move_vector):
        """
        function to move the scope by a certain distance vector, interpreted as a vector on the plane
        :param move_vector: tuple[float, float], using (x, y)
        :return: None
        """
        self._scope_position = (self._scope_position[0] + move_vector[0], self._scope_position[1] + move_vector[1])
        if self._scope_position[1] > 0:
            self._scope_position = (self._scope_position[0], 0)  # stick the coord to 0, if moved below x-axis
        self._set_border()

        self.redraw_plane(selector)

    def point_in_scope(self, point):
        """
        a function to figure out, if a given point is within the scope
        :param point: point on the plane
        :return: point, if in scope, False if not
        """
        # checking if point is confined in those edges
        if not (self._min_x_value <= point.real <= self._max_x_value):
            return False

        if not (self._min_y_value <= point.imag <= self._max_y_value):
            return False

        return point

    def get_crossing_of_arc_in_scope(self, centre_point, circle_radius):
        """
        calc the crossing points of scope and some arc
        :param centre_point: centre point of arc as complex
        :param circle_radius: radius of arc as float
        :return: list of 6 values, either a complex point or False:
                    (left_border_crossing, right_border_crossing, pos_down_border_crossing, neg_down_border_crossing,
                     pos_up_border_crossing, neg_up_border_crossing)
        """
        # init var
        min_x = centre_point.real - circle_radius
        max_x = centre_point.real + circle_radius
        min_y = 0
        max_y = circle_radius

        left_border_crossing = False
        right_border_crossing = False
        pos_down_border_crossing = False
        neg_down_border_crossing = False
        pos_up_border_crossing = False
        neg_up_border_crossing = False

        # borders:
        if self._min_x_value > min_x:
            left_border_crossing = complex(self._min_x_value,
                                           sqrt(sq(circle_radius) - sq(self._min_x_value - centre_point.real)))
        if max_x > self._max_x_value:
            right_border_crossing = complex(self._max_x_value,
                                            sqrt(sq(circle_radius) - sq(self._max_x_value - centre_point.real)))

        if self._min_y_value > min_y:
            pos_down_border_crossing = complex(sqrt(sq(circle_radius) - sq(self._min_y_value)) + centre_point.real,
                                               self._min_y_value)
            neg_down_border_crossing = complex(-sqrt(sq(circle_radius) - sq(self._min_y_value)) + centre_point.real,
                                               self._min_y_value)
        if max_y > self._max_y_value:
            pos_up_border_crossing = complex(sqrt(sq(circle_radius) - sq(self._max_y_value)) + centre_point.real,
                                             self._max_y_value)
            neg_up_border_crossing = complex(-sqrt(sq(circle_radius) - sq(self._max_y_value)) + centre_point.real,
                                             self._max_y_value)

        return (self.point_in_scope(neg_down_border_crossing), self.point_in_scope(left_border_crossing),
                self.point_in_scope(neg_up_border_crossing), self.point_in_scope(pos_up_border_crossing),
                self.point_in_scope(right_border_crossing), self.point_in_scope(pos_down_border_crossing))

    def get_crossing_of_line_in_scope(self, real_value):
        """
        returns the lower and upper crossing of the type 2 line with the scope of the plane
        :param real_value: the real value at which the type 2 line is located
        :return: tuple of (lower_crossing_point, upper_crossing_point)
        """
        return (complex(real_value, self._min_y_value), complex(real_value, self._max_y_value))

    def get_real_value_in_scope(self, real_value):
        return self._min_x_value < real_value < self._max_x_value


class BeltramiKleinDiskWidget:

    def __init__(self, x_coor, y_coor, size, color=(50, 50, 100)):
        self.position = (x_coor, y_coor)
        self._size = size  # tuple (x_size, y_size)
        self._radius = min(self._size[0], self._size[1])//2
        self._disk = BeltramiKleinDisk(self._radius, self.position[0] + self._radius, self.position[1] + self._radius, color)
        self.widget_identifying_color = color
        self.button_size = self._radius//4

        self.button_drag_drop = DragButton(self.position[0], self.position[1], self.button_size,
                                           img_blank, img_blank_pressed)
        self.button_rotate_counter_clock = UniversalSimpleButton(self.position[0] + self._size[0] - self.button_size,
                                                                 self.position[1], self.button_size,
                                                                 rotate_counter_clock,
                                                                 rotate_counter_clock, self._disk.change_orientation,
                                                                 False, QUARTER_PI/2)

        self.button_list = [self.button_drag_drop, self.button_rotate_counter_clock]

        bkdisk_list.append(self)

    def draw_border(self):
        """
        drawing a border around the widget for visual clarity
        :return: None
        """
        stroke(self.widget_identifying_color[0], self.widget_identifying_color[1], self.widget_identifying_color[2])
        fill(0, 0, 0)
        rectMode(CORNER)
        rect(self.position[0], self.position[1], self._size[0], self._size[1])

    def draw(self):
        """
        drawing the initial widget
        :return: None
        """
        self.draw_border()
        self._disk.draw()

    def redraw(self, selector=-1):
        """
        redrawing the disk of the widget
        :param selector: (int) if given redraws only 1 state, else redraws all states
        :return: None
        """
        self._disk.redraw(selector)

    def redraw_all(self, selector=-1):
        """
        redrawing the whole widget
        :param selector: (int) if given redraws only 1 state, else redraws all states
        :return: None
        """
        self.draw_border()
        self.redraw(selector)

        # redrawing buttons
        for button in self.button_list:
            button.update()

    def drag(self, vector):
        """
        changes the position of every element in the widget an of the widget itself by a vector
        :param vector: (int, int)
        :return: None
        """
        self.position = (self.position[0] + vector[0], self.position[1] + vector[1])
        self._disk.position = (self._disk.position[0] + vector[0], self._disk.position[1] + vector[1])
        for button in self.button_list:
            button.px += vector[0]
            button.py += vector[1]
        redraw_all()


class BeltramiKleinDisk:

    def __init__(self, radius, x_coor, y_coor, color):
        self._radius = radius
        self.position = (x_coor, y_coor)
        self._orientation = 0
        self._ident_color = color

    def draw_disk(self):
        """
        draws the empty disk
        :return: None
        """
        ellipseMode(CENTER)
        noStroke()
        fill(255, 255, 255)
        circle(self.position[0], self.position[1], 2*self._radius)

    def draw(self):
        """
        drawing the disk
        :return:
        """
        self.draw_disk()
        self._draw_orientation()

    def redraw(self, selector=-1):
        """
        drawing the disk with lines
        :return:
        """
        self.draw()
        if selector >= 0:
            if len(shapelist) > selector:  # check if selected shape exist for error catching
                for stored_line in shapelist[selector].lines:
                    stored_line.drawline_on_bk_disk(self)
            return None

        # draw lines on plane
        for shape in shapelist:
            for stored_line in shape.lines:
                stored_line.drawline_on_bk_disk(self)

    def _draw_orientation(self):
        """
        drawing the orientation onto the circle by coloring a line at the outside of the circle
        :return: None
        """
        stroke(self._ident_color[0], self._ident_color[1], self._ident_color[2])
        noFill()
        d = 20  # length of marker
        a = (origin[0] + math.cos(self._orientation)*(radius+5), origin[1] - math.sin(self._orientation)*(radius+5))
        b = (origin[0] + math.cos(self._orientation)*(radius+d), origin[1] - math.sin(self._orientation)*(radius+d))

        line(a[0], a[1], b[0], b[1])

    def get_point_absolute_pos(self, point):
        """
        convert relative point position to disk, into absolute point position on application
        :param point: complex point on disk
        :return: (int, int) of position in application
        """
        x_coor = point.real*self._radius + self.position[0]
        y_coor = point.imag*self._radius + self.position[1]
        return int(x_coor), int(y_coor)

    def set_orientation(self, value):
        """
        setter for orientation, if value outside [0, 2PI), it won't change the value
        :return: None
        """
        if 0 <= value < TWO_PI:
            self._orientation = value
            redraw_all(selector)

    def change_orientation(self, angle):
        """
        changes the orientation by an angle using the setter function
        :param angle: amount of change for angle
        :return: None
        """
        self.set_orientation((self._orientation + angle) % TWO_PI)

    def get_orientation(self):
        """
        getter for the _orientation param
        :return: _orientation in float
        """
        return self._orientation

    def calc_beltrami_klein_point(self, disk_point):
        """
        calculates the corresponding point on the beltrami-klein disk
        :param disk_point: (complex) point on the disk
        :return: corresponding point on the plane
        """
        disk_point = self.turn_circle(disk_point)
        temp = 1 + pow(disk_point.real, 2) + pow(disk_point.imag, 2)

        return complex(2*disk_point.real / temp, 2*disk_point.imag / temp)

    def turn_circle(self, disk_point):
        """
        interprets circle orientation on a given point
        :return: complex point
        """
        return cmath.exp(complex(0, self._orientation)) * disk_point


class Hyperbolic_Line:
    '''
    A class which represents a Line on the Pioncaré Disk and Plane
    '''

    def __init__(self, a, b, tcolor=[255, 0, 155], tpath=True):
        """
        input : a = (complex number) A Point on the Line
                b = (complex number) A point on the Line
                tcolor = (3 Element List of Int) represents RBG Code for Line color
                tpath = (boolean) will only draw the path from a to b
        Output: Instance of Hyperbolic_Line
        """
        self.disk_p1 = a  # point 1 of disk (used to create this line)
        self.disk_p2 = b  # point 2 of disk (used to create this line)

        # handling of choosing the same point twice
        if self.disk_p1 == self.disk_p2:
            self.disk_p2 = complex(0, 0)
            if self.disk_p1 == complex(0, 0):
                self.disk_p1 = complex(0.1, 0)

        self.plane_p1 = calc_plane_point(self.disk_p1, 0)
        self.plane_p2 = calc_plane_point(self.disk_p2, 0)

        # sorting plane_points using real value (x-coor), order is important for drawing in plane
        if self.plane_p1.real > self.plane_p2.real:
            self.plane_p1, self.plane_p2 = self.plane_p2, self.plane_p1
            self.disk_p1, self.disk_p2 = self.disk_p2, self.disk_p1

        self.plane_center = (0, 0)
        self.plane_center = self.plane_calc_center()
        self.color = tcolor
        self.typ_1_line = True  # marks if the line is of type 1 (Euclidean circle) or type 2 (Euclidean line) in plane
        self.disk_center = self.disk_calc_center()
        self.intersections = self.calc_circle_intersections()
        self.path = tpath

        # settings:
        self.plane_point_tracker = True
        self.disk_point_tracker = True
        redraw_all()

    def _reset_plane_data(self, plane):
        """
        setting up the values to draw on plane, same as in init,
        because it has to change based on plane orientation
        :return: None
        """
        self.plane_p1 = calc_plane_point(self.disk_p1, plane.get_orientation())
        self.plane_p2 = calc_plane_point(self.disk_p2, plane.get_orientation())

        self.typ_1_line = True
        # sorting plane_points using real value (x-coor), order is important for drawing in plane
        if self.plane_p1.real > self.plane_p2.real:
            self.plane_p1, self.plane_p2 = self.plane_p2, self.plane_p1
            self.disk_p1, self.disk_p2 = self.disk_p2, self.disk_p1

        self.plane_center = (0, 0)
        self.plane_center = self.plane_calc_center()

    def disk_calc_center(self):
        """
        calculates the center of the orto-cicle with given points
        """
        denominator = self.disk_p1 * self.disk_p2.conjugate() - self.disk_p1.conjugate() * self.disk_p2
        if denominator != complex(0, 0):
            numerator = (self.disk_p1 * (1 + abs(self.disk_p2) ** 2)
                         - self.disk_p2 * (1 + abs(self.disk_p1) ** 2))
            return numerator / denominator
        else:
            return None

    def plane_calc_center(self):
        """
        calculates the centre of the circle of the type 1 line, that sits on the x-axis
        :return: either None if the line is of type 2 and therefore not a circle,
                 or the centre point as a complex number.
        """
        if not abs(self.plane_p1.real - self.plane_p2.real) < 0.000001:  # if dif is lower than 10^-6 treat as if real values are the same
            # (p2_x^2 + p2_y^2 - p1_x^2 - p1_y^2) / (2(p2_x - p1_y))
            x = ((sq(self.plane_p2.real) + sq(self.plane_p2.imag) - sq(self.plane_p1.real) - sq(self.plane_p1.imag))
                 / (2 * (self.plane_p2.real - self.plane_p1.real)))  # see paper 14.1.24
            return complex(x, 0)
        else:
            self.typ_1_line = False
            return None

    def calc_circle_intersections(self):
        """
        calculates the intersections of the poincaré unit circle and the euclidian circle the hyperbolic line travels on
        (originaly calculated in drawline, but moved here for reuseability)
        :return: list of 2 complex numbers
        """
        if self.disk_center:
            # getting crossing point of the 2 circles
            radius = abs(self.disk_p1 - self.disk_center)
            return cicle_intersections(self.disk_center, radius, complex(0, 0), 1)
        else:
            # if it's a line crossing a circle, not 2 circles (special case when line crosses origin)
            # find line function
            if self.disk_p1.real - self.disk_p2.real == 0:
                return [complex(0, 1), complex(0, -1)]
            gradient = (self.disk_p1.imag - self.disk_p2.imag)/(self.disk_p1.real - self.disk_p2.real)
            y_axis_crossing = 0
            # solving with pq
            p = gradient*y_axis_crossing/(pow(gradient, 2)+1)
            q = (- 1)/(pow(gradient, 2) + 1)
            tmp = sqrt(pow(p, 2) - q)
            x_1 = p + tmp
            x_2 = p - tmp
            return [complex(x_1, gradient*x_1), complex(x_2, gradient*x_2)]

    def drawline(self):
        """
        draws the hyperbolic line on the disk
        :return: None
        """
        # draw line on disk
        if not self.disk_center or abs(self.disk_center) > 2**45:  # special case where the line is an euclidian line
            if self.path:
                # get absolute points in app
                a = interpret_point(self.disk_p1)
                b = interpret_point(self.disk_p2)
            else:
                a, b = interpret_point(self.intersections[0]), interpret_point(self.intersections[1])

            # draw line
            stroke(self.color[0], self.color[1], self.color[2])
            line(a[0], a[1], b[0], b[1])

        # draw circle on disk
        else:
            global radius
            r = abs(self.disk_p1 - self.disk_center)
            p = interpret_point(self.disk_center)
            if self.path == True:
                radians1 = cmath.phase(self.disk_p1 - self.disk_center)
                radians2 = cmath.phase(self.disk_p2 - self.disk_center)
            else:
                radians1 = cmath.phase(self.intersections[0] - self.disk_center)
                radians2 = cmath.phase(self.intersections[1] - self.disk_center)
            # check if we got the smaller arc
            if abs(radians1 - radians2) > PI:
                if radians1 < radians2:
                    radians1 += 2 * PI
                else:
                    radians2 += 2 * PI
            if radians1 > radians2:
                radians1, radians2 = radians2, radians1
            stroke(self.color[0], self.color[1], self.color[2])
            noFill()
            arc(p[0], p[1], r * 2 * radius, r * 2 * radius, radians1, radians2)

        if self.disk_point_tracker:
            rectMode(CENTER)
            plane_p1_tmp = interpret_point(self.disk_p1)
            plane_p2_tmp = interpret_point(self.disk_p2)
            rect(plane_p1_tmp[0], plane_p1_tmp[1], 4, 4)
            rect(plane_p2_tmp[0], plane_p2_tmp[1], 4, 4)

    def drawline_on_plane(self, plane):
        """
        draws the line on a given plane object
        :param plane: HyperbolicPlane, on which to draw the line
        :return: None
        """

        # draw settings
        stroke(self.color[0], self.color[1], self.color[2])
        noFill()
        self._reset_plane_data(plane)

        if self.typ_1_line:
            border_points = []  # the points that signal the start and end of where to draw

            # check if path or line
            if self.path:
                if plane.point_in_scope(self.plane_p1):
                    border_points.append(self.plane_p1)
                if plane.point_in_scope(self.plane_p2):
                    border_points.append(self.plane_p2)
            else:
                # calc crossing points with x-axis (center - radius and center + radius) and treat it as start/end point
                if plane.point_in_scope(complex(self.plane_center.real - abs(self.plane_p1 - self.plane_center), 0)):
                    border_points.append(complex(self.plane_center.real - abs(self.plane_p1 - self.plane_center), 0))
                if plane.point_in_scope(complex(self.plane_center.real + abs(self.plane_p1 - self.plane_center), 0)):
                    border_points.append(complex(self.plane_center.real + abs(self.plane_p1 - self.plane_center), 0))

            # calc crossing points
            crossing_result = plane.get_crossing_of_arc_in_scope(self.plane_center,
                                                                 abs(self.plane_p1 - self.plane_center))

            for point in crossing_result:
                if point:
                    # adding points if crossing point exists
                    if plane.point_in_scope(point) and (not self.path or self.plane_p1.real < point.real < self.plane_p2.real):
                        # adding points if in scope and (line, or on path)
                        border_points.append(point)

            border_points = sorted(border_points, key=return_real_value)
            """# tmp drawing crossing points todo: remove
            for crossing in crossing_result:
                if crossing:
                    rectMode(CENTER)
                    cross_tmp = plane.get_point_relative_to_plane(crossing)
                    rect(cross_tmp[0], cross_tmp[1], 4, 4)"""

            line_list = []
            last_crossing = False
            # creating a list of lines part of the hyperbolic line, that should be visible
            for crossing in border_points:
                if last_crossing and crossing:
                    line_list.append((last_crossing, crossing))
                    last_crossing = False
                    continue
                if crossing:
                    last_crossing = crossing

            # getting centre of circle point
            m_point = plane.get_point_relative_to_plane(self.plane_center)

            # calc diameter of the circle and fit to scope
            distance_vector = plane.get_scope_relative_scaling(self.plane_p1 - self.plane_center)
            plane_radius = abs(complex(distance_vector[0], distance_vector[1]))
            plane_diameter = (plane_radius * 2, plane_radius * 2)

            for lines in line_list:
                # calc radians
                plane_radian = (-cmath.phase(lines[1] - self.plane_center),
                                -cmath.phase(lines[0] - self.plane_center))  # gonna be needed for parts of the semi circle (lines)
                if plane_radian[0] > plane_radian[1]:
                    plane_radian = (plane_radian[1], plane_radian[0])  # swap pos, so the first is the smaller number

                # draw arc
                arc(m_point[0], m_point[1], plane_diameter[0], plane_diameter[1], plane_radian[0], plane_radian[1])

        elif plane.get_real_value_in_scope(self.plane_p1.real):  # type 2 line
            a = self.plane_p1
            b = self.plane_p2
            # get border crossing points:
            lower_bound, upper_bound = plane.get_crossing_of_line_in_scope(self.plane_p1.real)

            # if a point is out of bounds, set it to the end of scope
            if lower_bound.imag > a.imag:
                a = lower_bound
            if upper_bound.imag < a.imag:
                a = upper_bound
            if lower_bound.imag > b.imag:
                b = lower_bound
            if upper_bound.imag < b.imag:
                b = upper_bound

            # abs point calc
            if not self.path:
                a = lower_bound
                b = upper_bound


            a = plane.get_point_relative_to_plane(a)
            b = plane.get_point_relative_to_plane(b)
            if not a == b:  # draw a line, exept if points are both out of scope on the same side
                line(a[0], a[1], b[0], b[1])

        # draw points, used to create the line on plane
        if self.plane_point_tracker:
            rectMode(CENTER)
            if plane.point_in_scope(self.plane_p1):
                plane_p1_tmp = plane.get_point_relative_to_plane(self.plane_p1)
                rect(plane_p1_tmp[0], plane_p1_tmp[1], 4, 4)
            if plane.point_in_scope(self.plane_p2):
                plane_p2_tmp = plane.get_point_relative_to_plane(self.plane_p2)
                rect(plane_p2_tmp[0], plane_p2_tmp[1], 4, 4)

    def drawline_on_bk_disk(self, disk):
        """
        function to draw line on a beltrami-klein disk
        :param disk: an BeltramiKleinDisk
        :return: None
        """
        bk_point_1 = disk.calc_beltrami_klein_point(self.disk_p1)
        bk_point_2 = disk.calc_beltrami_klein_point(self.disk_p2)
        # draw settings
        stroke(self.color[0], self.color[1], self.color[2])
        noFill()
        if self.path:
            a = bk_point_1
            b = bk_point_2
        else:
            # for lines, we have the same intersections, as with the poincaré disk,
            #  meaning we can reuse the calculations from there
            a = disk.turn_circle(self.intersections[0])
            b = disk.turn_circle(self.intersections[1])
        a = disk.get_point_absolute_pos(a)
        b = disk.get_point_absolute_pos(b)
        line(a[0], a[1], b[0], b[1])

        # draw points, used to create the line on plane
        if self.plane_point_tracker:
            abs_p1 = disk.get_point_absolute_pos(bk_point_1)
            abs_p2 = disk.get_point_absolute_pos(bk_point_2)
            rectMode(CENTER)
            rect(abs_p1[0], abs_p1[1], 4, 4)
            rect(abs_p2[0], abs_p2[1], 4, 4)

    def mirror(self, z):
        """
        input : z = (complex) point to be mirrored
        output:  (complex) mirrored point z'
        """
        if self.disk_center:
            # mirroring on a circle
            numerator = self.disk_center * z.conjugate() - 1
            denominator = z.conjugate() - self.disk_center.conjugate()
            return numerator / denominator
        else:
            # mirroring on a line
            v = self.disk_p1 - self.disk_p2
            projection = (z.real*v.real + z.imag*v.imag)*v/pow(abs(v), 2)
            return 2*projection - z

    def translate(self, a, tcolor=[255, 0, 0]):
        """
        input : a = (complex) vector to translate line
        tcolor = (3 Element List of Int) represents RBG Code for Line color the translated line
        output: (instance of Hyperbolic_Line) translated line with vector a
        """
        return Hyperbolic_Line(P_disc.translation(self.disk_p1, a), P_disc.translation(self.disk_p2, a), tcolor,
                               self.path)

    def intuitive_disk_translate(self, translation_function):
        """
        this function should get a function from the P_disk.intuitive_translation_constructor and creates a new
        hyperbolic line object that has been translated
        :param translation_function: a function that translates points from the disk onto other points on the disk
        :return: Hyperbolic_Line
        """
        return Hyperbolic_Line(translation_function(self.disk_p1), translation_function(self.disk_p2))

    """    
    # function zum testen der genauigkeit der abbildung
    def test(self):
        num_list = [complex(0.5, 0.5), complex(0.7, 0.7), complex(0.9, 0.9), complex(0.1, 0.1), complex(0.1, 0.9),
                    complex(0.3, 0.7), complex(0.8, 0.5), complex(0.4, 0.9), complex(0.001, 0.999), complex(0.99999, 0.00001)]

        for x in num_list:
            y = self.calc_plane_point(x)
            print(x, " -> ", y, "\n")
        for x in num_list:
            y = self.calc_plane_point(-x)
            print(-x, " -> ", y, "\n")
    """


def calc_plane_point(disk_point, angle):
    """
    uses the Cayley-transformation to calculate a point on the plane, based on a point on the disk
    :param disk_point: (complex) point on the disk
    :param angle: the orientation of the scope, values of [0, 2PI) are allowed
    :return: corresponding point on the plane
    """
    disk_point = cmath.exp(complex(0, angle)) * disk_point  # e^theta * z
    return complex(0, 1) * (1 + disk_point) / (1 - disk_point)  # i · (1 + z) / (1 − z)


def return_real_value(value):
    """
    function used as a key to sort numbers by their real value, if possible use value.real everyother time
    :param value: a complex number
    :return: a float amounting to the real value of the complex number
    """
    return value.real


class Hyperbolic_Shape:
    '''
    A class filled with instances of Hyberbolic_line to group them into one Shape
    '''

    def __init__(self, *input):
        """
        input : any Number of instances of Hyberbolic_line
        Output: Instance of Hyperbolic_Shape
        """
        # self.n = 0
        # self.edges = []
        self.lines = []
        shapelist.append(self)
        for l in input:
            self.include(l)

    def include(self, l):
        """
        Includes a Hyberbolic_line instance into the Shape
        input : l = Hyberbolic_line instance
        """
        self.lines.append(l)
        # self.edges.extend([l.p1,l.p2])
        # self.n += 1

    def tesselation(self, m=10, ignore_first=False):
        """
        tesselation of the Pioncare Disk with given Shape
        input : m = (int) Loopsteps left until the recursion stops
                ignore_first = (boolean) should be False if called , used to skip the first line of the used shape to stop remirroring already calcaulted and drawn shapes
        """
        if m > 0:
            for l in self.lines:
                if ignore_first == True:
                    ignore_first = False
                    continue
                mirrored = Hyperbolic_Shape(l)
                for j in self.lines:
                    if l != j:
                        mirrored.include(Hyperbolic_Line(l.mirror(j.disk_p1), l.mirror(j.disk_p2), j.color, j.path))
                        mirrored.tesselation(m - 1, True)


class Button:
    """
    represent a Menu Button, should only be inherited from, has no action on press
    """

    def __init__(self, x, y, tpic, tpressedpic):
        """
        input : x = (int) x-Position of the Button on a 20x20 grid
                y = (int) y-Position of the Button on a 20x20 grid
                tpic = (processing pic) pictures of the unpressed button allready loaded in load_images()
                tpressedpic = (processing pic) pictures of the pressed button allready loaded in load_images()
        Output: Instance of Button
        """
        global buttonlist
        self.px = float(x)
        self.py = float(y)
        self.pic = tpic
        self.pressed_pic = tpressedpic
        self.pressed = False
        buttonlist.append(self)

        if width//20 > height//10:
            self.b_size = height//10
        else:
            self.b_size = width//20

        while self.px < 0:
            self.px += width//self.b_size
        while self.py < 0:
            self.py += height//self.b_size

    def update(self):
        """
        method to draw the button with the right image
        """
        # introducing b_size, so that the images don't get stretched the same way, the window does
        # todo: for some reason I cant define b_size already in the init -> reason found, since super() doesnt work, there is no b_size in the inherited classes
        if self.pressed:
            image(self.pressed_pic, self.b_size * self.px, self.b_size * self.py, self.b_size,
                  self.b_size)
        else:
            image(self.pic, self.b_size * self.px, self.b_size * self.py, self.b_size,
                  self.b_size)

    def action(self):
        """
        method to handle the button press
        """
        return "Nothing todo here"

    def check(self, pos):
        """
        method to check if a point(mouseclick) is inside this button, so the button can be considered pressed
        input : pos ( 2 Element List of Integer) Position of the Point on the Canvas
        """
        if (self.px * self.b_size <= pos[0] < (self.px + 1) * self.b_size and self.py *
                self.b_size <= pos[1] < (self.py + 1) * self.b_size):
            self.action()
            return True  # indicator to break loops


class UniversalSimpleButton(Button):
    def __init__(self, x, y, b_size, tpic, tpressedpic, action_function, switch, *args):
        """
        :param x: x_coordinate in absolute value
        :param y: y_coordinate in absolute value
        :param tpic: non pressed button image
        :param tpressedpic: pressed button image
        :param action_function: a function without param that should be run when the button is clicked
        :param switch: bool, if the button should switch between pressed and not pressed state with each click
        :param args: arguments to give the button function, not useable on switch
        """
        self.px = x
        self.py = y
        self.pic = tpic
        self.pressed_pic = tpressedpic
        self.pressed = False
        self.b_size = b_size

        self.switch = switch
        self._function = action_function
        self._arguments = args[0]
        self.update()

    def action(self):
        if self.switch:
            self._function()
        else:
            self._function(self._arguments)
        self.update()

    def update(self):
        if self.pressed:
            image(self.pressed_pic, self.px, self.py, self.b_size, self.b_size)
        else:
            image(self.pic, self.px, self.py, self.b_size, self.b_size)

    def check(self, pos):
        """
        method to check if a point(mouseclick) is inside this button, so the button can be considered pressed
        input : pos ( 2 Element List of Integer) Position of the Point on the Canvas
        """
        if self.px <= pos[0] < self.px + self.b_size and self.py <= pos[1] < self.py + self.b_size:
            if self.switch:
                self.pressed = not self.pressed
            self.action()
            return True  # indicator to break loops


class Number_button(Button):
    """
    represent a Number Button
    """

    def __init__(self, x, y, tpic, tpressedpic, tnumber):
        """
        for whatever reason super() seems not to work on this jyhton interpreter , so here we go again
        input : same as in Button
                tumber (int) Number associated with this button
        Output: Instance of Number_button
        """
        self.px = x
        self.py = y
        self.pic = tpic
        self.pressed_pic = tpressedpic
        self.pressed = False
        self.number = tnumber

        if width // 20 > height // 10:
            self.b_size = height // 10
        else:
            self.b_size = width // 20
        while self.px < 0:
            self.px += width//self.b_size
        while self.py < 0:
            self.py += height//self.b_size

        buttonlist.append(self)

    def action(self):
        """
        method to handle the button press , handles animations and calls number_pressed()
        """
        global lastnumberpressed
        lastnumberpressed.pressed = False
        lastnumberpressed.update()
        lastnumberpressed = self
        self.pressed = True
        self.update()
        number_pressed(self.number)


class Tool_button(Button):

    def __init__(self, x, y, tpic, tpressedpic, ttool):
        """
        for whatever reason super() seems not to work on this jyhton interpreter , so here we go again
        input : same as in Button
                ttool (int) Number of the Tool associated with this button
        Output: Instance of Tool_button
        """
        self.px = x
        self.py = y
        self.pic = tpic
        self.pressed_pic = tpressedpic
        self.toolnumber = ttool
        self.pressed = False

        if width // 20 > height // 10:
            self.b_size = height // 10
        else:
            self.b_size = width // 20
        while self.px < 0:
            self.px += width//self.b_size
        while self.py < 0:
            self.py += height//self.b_size

        buttonlist.append(self)

    def action(self):
        """
        method to handle the button press, handles animations,
        checks if the unit cicle need to be redrawn and calls change_mode()
        """
        global lasttoolpressed
        global redraw_lines
        clear_textbox()
        if redraw_lines == True:
            redraw_all()
        lasttoolpressed.pressed = False
        lasttoolpressed.update()
        self.pressed = True
        self.update()
        lasttoolpressed = self
        change_mode(self.toolnumber)


class DragButton(Button):
    """
    A button specificly designed to drag'n'drop a widget
    """

    def __init__(self, x, y, bsize, tpic, tpressedpic):
        """
        input : x = (int) absolute x position
                y = (int) absolute y position
                bsize = (int) size of the img in pixel
                tpic = (processing pic) pictures of the unpressed button already loaded in load_images()
                tpressedpic = (processing pic) pictures of the pressed button already loaded in load_images()
        Output: Instance of Button
        """
        global buttonlist
        self.px = float(x)
        self.py = float(y)
        self.pic = tpic
        self.pressed_pic = tpressedpic
        self.pressed = False  # True while holding mouse button down
        self.b_size = bsize

        self.stored_mouse_position = (0, 0)  # mouse position when clicked

    def action(self):
        self.drag()
        self.update()

    def update(self):
        if self.pressed:
            image(self.pressed_pic, self.px, self.py, self.b_size, self.b_size)
        else:
            image(self.pic, self.px, self.py, self.b_size, self.b_size)

    def drag(self):
        """
        saves mouse position on click
        :return: None
        """
        self.pressed = True
        self.stored_mouse_position = (mouseX, mouseY)
        return None

    def drop(self):
        """
        gets called, by mouseReleased(), returns a delta value to calc new position of widget,
        does nothing when button wasn't pressed
        :return: delta of mouse position (int, int)
        """
        if self.pressed:
            self.pressed = False
            self.update()
            return mouseX - self.stored_mouse_position[0], mouseY - self.stored_mouse_position[1]
        return 0, 0

    def check(self, pos):
        """
        method to check if a point(mouseclick) is inside this button, so the button can be considered pressed
        input : pos ( 2 Element List of Integer) Position of the Point on the Canvas
        """
        if self.px <= pos[0] < self.px + self.b_size and self.py <= pos[1] < self.py + self.b_size:
            self.action()
            return True  # indicator to break loops


# do not change Variables
origin = [screen_x // 4, screen_y // 2]  # origin of the unit cicle
radius = float(min(screen_x, screen_y)) / 3  # radius of the unit cicle
stored = complex(0, 0)  # stored the last clicked point on the unit circle
clicks = 0  # hold the amounts of clicks
mode = 1  # stores the menu mode
shapelist = []  # list of all hyberbolic shapes
buttonlist = []  # list of all buttons
lastnumberpressed = None  # keeps the last number pressed for animation
lasttoolpressed = None  # keeps the last tool pressed for animation
selector = -1  # holds which layer is selected
to_be_mirrored = None  # holds the line which has to be mirrored
amount = None  # holds the amount of interations
redraw_lines = False  # keeps track if all lines have to be redrawn on tool change
typing = ""  # holds the keyboard inputs
input = "input"  # holds the keyboard input after ENTER

plane_list = []  # a collection of every HPlaneWidget
bkdisk_list = []  # a collection of every BeltramiKleinDiskWidget


def interpret_point(input):
    global radius
    global origin
    '''
    Interpret a complex number z into a touple of R^2 , does a basis change and translate
    its location based on the cicle midpoint "origin" and radius of the cicle
    
    For example if midpoint M(m1,m2) with  radius r, then z = x+iy will be
     mapped to [x*r + m1, y*r + m2]
     
    input: input : complex number or List(2 Elements of floats)
    output: P: List of 2 real numbers forming a point in R^2 or a complex number
    '''
    if type(input) == complex:  # from C to R^2
        return [input.real * radius + origin[0], input.imag * radius + origin[1]]

    elif type(input) == list:  # from R^2 to C
        return complex((input[0] - origin[0]) / radius, (input[1] - origin[1]) / radius)

    else:
        print("invalid Type")
        return None


def cicle_intersections(m1, r1, m2, r2, d=None):
    """
    function to calculate the intersections of two circles, optional distance (d) to cut calculations
    m1 = (complex) Midpoint of circle 1
    r1 = (float) radius of circle 1
    m2 = (complex) Midpoint of circle 2
    r2 = (float) radius of circle 2
    d = (float) parameter for iteral calls with same distance to cut calculations
    return: (list of complex) both intersection of circle 1 and 2.
    """
    if d == None:
        d = abs(m1 - m2)
    if d <= r1 + r2:
        v = (m2 - m1) / d  # normalized direction vector M1M2
        nv = complex(- v.imag, v.real)
        # negative reciprocal of v to get its normal vector
        if r1 == r2:  # common case, optimated calculation
            t = d / 2  # distance from M1 to intersection of M1M2 and P1P2
        else:
            t = (r1 ** 2 - r2 ** 2 + d ** 2) / (2 * d)
        l = sqrt(r1 ** 2 - t ** 2)  # length of the path from P1 to P2
        s = m1 + t * v  # move to intersection from M1M2 and P1P2
        return [s + l * nv, s - l * nv]  # move up to P1, move down to P2
    else:  # sum of radii to small
        return None


def set_selector(number):
    """
    setting the selector
    :param number: (int) number to set selector to
    :return: None
    """
    global buttonlist
    global selector
    for button in buttonlist:
        if isinstance(button, Number_button) and button.number == number:
            button.action()
            return
    selector = int(number)


# random line generation
import random


def random_two_points_on_disk():
    """
    Random number generator on unitcircle
    :return: tuple of 2 complex numbers
    """
    t, u = random.random(), random.random()
    point_1 = complex(math.sqrt(t) * math.cos(2 * math.pi * u), math.sqrt(t) * math.sin(2 * math.pi * u))
    t, u = random.random(), random.random()
    point_2 = complex(math.sqrt(t) * math.cos(2 * math.pi * u), math.sqrt(t) * math.sin(2 * math.pi * u))

    return point_1, point_2


def generate_n_random_lines(n=10):
    """
    gerneates n lines, and appends it to random selectors
    :param n: (int) number of iterations, every iteration creates 1 line
    :return: None
    """
    for _ in range(n):
        selec = random.randrange(0, len(shapelist))
        a, b = random_two_points_on_disk()
        shapelist[selec].include(Hyperbolic_Line(a, b, tcolor=[0, 200, 0], tpath=bool(random.getrandbits(1))))



# Hud functions

def load_images():
    """
    loads in all icon images , this can only happend inside setup() and still needs to be declared to be global so they can be seen in draw().
    else the image have to be loaded each time they should be drawn.
    """
    # ui button pics
    global line_pic
    global line_pic_pressed
    global path_pic
    global path_pic_pressed
    global mirror_pic
    global mirror_pic_pressed
    global transpose_pic
    global transpose_pic_pressed
    global equal_pic
    global equal_pic_pressed
    global tesselate_pic
    global tesselate_pic_pressed
    global select_pic
    global select_pic_pressed
    global screen_pic
    global screen_pic_pressed
    global trash_pic
    global trash_pic_pressed
    global zero_pic
    global zero_pic_pressed
    global one_pic
    global one_pic_pressed
    global two_pic
    global two_pic_pressed
    global three_pic
    global three_pic_pressed
    global four_pic
    global four_pic_pressed
    global five_pic
    global five_pic_pressed
    global six_pic
    global six_pic_pressed
    global seven_pic
    global seven_pic_pressed
    global eight_pic
    global eight_pic_pressed
    global nine_pic
    global nine_pic_pressed

    # hyperbolic plane widget button pics
    global plus_image
    global minus_image
    global up_image
    global down_image
    global left_image
    global right_image
    global plus_pressed_image
    global minus_pressed_image
    global up_pressed_image
    global down_pressed_image
    global left_pressed_image
    global right_pressed_image
    global rotate_counter_clock
    global img_axis
    global img_axis_pressed
    global img_coordinate
    global img_coordinate_pressed
    global img_grid
    global img_grid_pressed
    global img_blank
    global img_blank_pressed
    global img_none
    global img_none_pressed


    line_pic = loadImage("line.png")
    line_pic_pressed = loadImage("line_pressed.png")
    path_pic = loadImage("path.png")
    path_pic_pressed = loadImage("path_pressed.png")
    mirror_pic = loadImage("mirror.png")
    mirror_pic_pressed = loadImage("mirror_pressed.png")
    transpose_pic = loadImage("transpose.png")
    transpose_pic_pressed = loadImage("transpose_pressed.png")
    equal_pic = loadImage("equal.png")
    equal_pic_pressed = loadImage("equal_pressed.png")
    tesselate_pic = loadImage("tesselate.png")
    tesselate_pic_pressed = loadImage("tesselate_pressed.png")
    select_pic = loadImage("select.png")
    select_pic_pressed = loadImage("select_pressed.png")
    screen_pic = loadImage("screen.png")
    screen_pic_pressed = loadImage("screen_pressed.png")
    trash_pic = loadImage("trash.png")
    trash_pic_pressed = loadImage("trash_pressed.png")
    zero_pic = loadImage("0.png")
    zero_pic_pressed = loadImage("0_pressed.png")
    one_pic = loadImage("1.png")
    one_pic_pressed = loadImage("1_pressed.png")
    two_pic = loadImage("2.png")
    two_pic_pressed = loadImage("2_pressed.png")
    three_pic = loadImage("3.png")
    three_pic_pressed = loadImage("3_pressed.png")
    four_pic = loadImage("4.png")
    four_pic_pressed = loadImage("4_pressed.png")
    five_pic = loadImage("5.png")
    five_pic_pressed = loadImage("5_pressed.png")
    six_pic = loadImage("6.png")
    six_pic_pressed = loadImage("6_pressed.png")
    seven_pic = loadImage("7.png")
    seven_pic_pressed = loadImage("7_pressed.png")
    eight_pic = loadImage("8.png")
    eight_pic_pressed = loadImage("8_pressed.png")
    nine_pic = loadImage("9.png")
    nine_pic_pressed = loadImage("9_pressed.png")

    plus_image = loadImage("plus.png")
    minus_image = loadImage("minus.png")
    up_image = loadImage("up.png")
    down_image = loadImage("down.png")
    left_image = loadImage("left.png")
    right_image = loadImage("right.png")
    plus_pressed_image = loadImage("plus_pressed.png")
    minus_pressed_image = loadImage("minus_pressed.png")
    up_pressed_image = loadImage("up_pressed.png")
    down_pressed_image = loadImage("down_pressed.png")
    left_pressed_image = loadImage("left_pressed.png")
    right_pressed_image = loadImage("right_pressed.png")
    rotate_counter_clock = loadImage("Rotate_Counterclockwise.png")
    img_axis = loadImage("axis.png")
    img_axis_pressed = loadImage("axis_pressed.png")
    img_coordinate = loadImage("coordinate.png")
    img_coordinate_pressed = loadImage("coordinate_pressed.png")
    img_grid = loadImage("grid.png")
    img_grid_pressed = loadImage("grid_pressed.png")
    img_blank = loadImage("blank.png")
    img_blank_pressed = loadImage("blank_pressed.png")
    img_none = loadImage("none.png")
    img_none_pressed = loadImage("none_pressed.png")


def setup_hud():
    """
    function to create all Buttons instances for the user Interface
    """
    global lastnumberpressed
    global lasttoolpressed
    B_line = Tool_button(0, 0, line_pic, line_pic_pressed, 0)
    B_path = Tool_button(1, 0, path_pic, path_pic_pressed, 1)
    B_mirror = Tool_button(2, 0, mirror_pic, mirror_pic_pressed, 2)
    B_transpose = Tool_button(3, 0, transpose_pic, transpose_pic_pressed, 3)
    B_equal = Tool_button(1, 1, equal_pic, equal_pic_pressed, 4)
    B_tesselate = Tool_button(2, 1, tesselate_pic, tesselate_pic_pressed, 5)
    B_select = Tool_button(0, 2, select_pic, select_pic_pressed, 6)
    B_screen = Tool_button(1, 2, screen_pic, screen_pic_pressed, 7)
    B_trash = Tool_button(0, 3, trash_pic, trash_pic_pressed, 8)
    # Path Tool is selected on Programm start, fix animation
    lasttoolpressed = B_path
    lasttoolpressed.pressed = True

    B_None = Number_button(-5, 0, img_none, img_none_pressed, -1)
    B_0 = Number_button(-4, 0, zero_pic, zero_pic_pressed, 0)
    B_1 = Number_button(-3, 0, one_pic, one_pic_pressed, 1)
    B_2 = Number_button(-2, 0, two_pic, two_pic_pressed, 2)
    B_3 = Number_button(-1, 0, three_pic, three_pic_pressed, 3)
    B_4 = Number_button(-3, 1, four_pic, four_pic_pressed, 4)
    B_5 = Number_button(-2, 1, five_pic, five_pic_pressed, 5)
    B_6 = Number_button(-1, 1, six_pic, six_pic_pressed, 6)
    B_7 = Number_button(-2, 2, seven_pic, seven_pic_pressed, 7)
    B_8 = Number_button(-1, 2, eight_pic, eight_pic_pressed, 8)
    B_9 = Number_button(-1, 3, nine_pic, nine_pic_pressed, 9)
    # TODO: Layer 0 is selected on Programm start, fix animation
    lastnumberpressed = B_0
    lastnumberpressed.pressed = True


def update_hud():
    """
    function to draw all buttons in buttonlist
    """
    for button in buttonlist:
        button.update()


def clear_textbox():
    """
    function to clear the textbox
    """
    stroke(255, 255, 255)
    fill(255, 255, 255)
    rect(0.5 * width / 20, 18.5 * height / 20, 3.5 * width / 20, height / 20, 10)


def redraw_all(select=None):
    """
    function to redraw all shapes in shapelist, as well as all widgets and the unit circle.
    Used for the select , photo and delete mode and many other redraw functions.
    """
    global shapelist
    global raidus
    global origin
    global selector
    if not isinstance(select, int):
        select = selector
    background(0, 0, 0)
    # clear poincare Unit Cicle
    fill(255)
    stroke(255, 255, 255)
    circle(origin[0], origin[1], 2*radius+50)
    noStroke()
    fill(0, 0, 0)
    circle(origin[0], origin[1], 2*radius+45)
    fill(255)
    stroke(255, 255, 255)
    circle(origin[0], origin[1], 2 * radius)

    # draw poincare circle lines
    if select >= 0:
        while len(shapelist) <= select:
            newshape = Hyperbolic_Shape()
        for l in shapelist[select].lines:
            l.drawline()
    else:
        for s in shapelist:
            for l in s.lines:
                l.drawline()
    # draw bkdisks + lines
    for disk in bkdisk_list:
        disk.redraw_all(select)
    # draw poincare plane + plane lines
    for hplane in plane_list:
        hplane.redraw_all(select)
    for button in buttonlist:
        button.update()
    redraw_lines = False  # set redraw flag off


def interpret_input(input):
    """
    function to interpret the input given. Only numbers are actually used , but takes str for implementation for commands of some sort
    input: input (str) to check
    """
    global selector
    global mode
    is_int = True
    is_float = True
    # check if the string can be converted to int and float
    try:
        int(input)
    except ValueError:
        is_int = False
    try:
        float(input)
    except ValueError:
        is_float = False

    if is_float == True and mode == 9:  # second instance of equal shape, input via keyboard
        P_disc.equilateral_Shape(selector, abs(float(input)))
    elif is_int == True and mode != 9:
        if int(input) > 9:  # number input higher the any Number Button
            lastnumberpressed.pressed = False
            lastnumberpressed.update()
            number_pressed(int(input))
        else:
            for button in buttonlist:  # check which button corresponds to given number for animation
                if isinstance(button, Number_button) == True and button.number == int(input):
                    button.action()
                    break


def change_mode(n):
    """
    function to react to mode changes , mostly used to give user feedback , but also used to get 1 or 2 click actions
    input: n (int) the number mode to replay to
    """
    global mode
    global shapelist
    global selector
    global buttonlist
    # Check for double presses
    if mode == 8:  # delete button pressed twice: eradicte everything!
        del shapelist
        shapelist = []
        mode = 1
        for button in buttonlist:
            if isinstance(button, Number_button) and button.number == -1:
                button.action()
                break
        new_shape = Hyperbolic_Shape()

    mode = n
    if mode == 2:  # Mirror case
        textbox_print("Mirror which layer?")
    elif mode == 3:  # transaltion case
        textbox_print("Move which layer?")

    elif mode == 4:  # first instance equal shape case
        textbox_print("How many edges?")

    elif mode == 5:  # tesselation case
        textbox_print("How many iterations?")

    elif mode == 6:  # selecting case
        textbox_print("Show which layer?")

    elif mode == 7:  # screenshot without user interface
        # can't overdraw button wuhtout overdrawing the unit circle , so clear the whole canvas
        background(0)
        fill(255)
        stroke(255, 255, 255)
        circle(origin[0], origin[1], 2 * radius)
        redraw_all()
        id_string = (str(year()) + "-" + str(month()) + "-" + str(day()) + "_"
                     + str(hour()) + "-" + str(minute()) + "-" + str(second()))
        save("Screenshot_" + id_string + ".png")
        update_hud()
        textbox_print("Photo taken!")

    elif mode == 8:  # delete case
        textbox_print("Delete which layer?")


def number_pressed(number):
    """
    function to handle number input for modes which uses the number inputs for everything else then layer selection
    input: number (int) the number to act on
    """
    global to_be_mirrored
    global selector
    global mode
    global amount
    selector = number
    # Statecheck , which mode is active?
    if mode == 2:  # mirror mode
        if to_be_mirrored == None and len(shapelist) >= selector:  # user want this layer to me mirrored
            to_be_mirrored = selector  # set layer to be mirrored
            textbox_print("On which layer?")
        elif len(shapelist) >= selector:  # user want this layer to be the mirror
            newshape = Hyperbolic_Shape()
            for l in shapelist[to_be_mirrored].lines:
                for j in shapelist[selector].lines:
                    newshape.include(Hyperbolic_Line(j.mirror(l.p1), j.mirror(l.p2), l.color, l.path))
            to_be_mirrored = None
            shapelist.append(newshape)
        else:
            textbox_print("Layer is Empty!")
    elif mode == 4:  # equal shape mode, n has been selected
        if selector > 2:
            mode = 9  # set to second instance of equal shape
            textbox_print("Radius?")
        else:  # errorcatch n<3
            textbox_print("Needs at least 3!")
    elif mode == 5:  # tesselation Case
        if amount == None:  # get n
            amount = selector
            textbox_print("Which layer?")
        elif len(shapelist) >= selector:
            shapelist[selector].tesselation(amount)
            amount == None
    elif mode == 6:  # selecting Case
        selector = number
        redraw_all(selector)
    elif mode == 8 and selector != -1:  # Delete Case
        if len(shapelist) > selector:  # check if selected shape exist for error catching
            del shapelist[selector].lines
            shapelist[selector].lines = []
            redraw_all()
    else:  # selecting layer (same as case mode 6, as that one is legacy)
        redraw_all(selector)
        

def textbox_print(string):
    """
    function to print strings in the textbox
    input: string (string) message to be printed in the textbox
    """
    clear_textbox()
    fill(0, 0, 0)  # change fontcolor to black
    text(string, 0.6 * width / 20, 18.8 * height / 20, 3.5 * width / 20, height / 20)


# Processing Functions        

def setup():
    """
    Runs on once on programm start, processing function for setup can only be run here
    """
    global screen_x
    global screen_y
    global origin
    global radius
    size(screen_x, screen_y)
    textSize(screen_x / 55)
    background(0)
    fill(255)
    stroke(255, 255, 255)
    circle(origin[0], origin[1], 2 * radius)
    load_images()
    setup_hud()
    update_hud()
    clear_textbox()

    hplane_widget_one = HyperbolicPlaneWidget(x_coor=screen_x//2, y_coor=20, size=(400, 260), widget_color=(200, 200, 0))
    #hplane_widget_two = HyperbolicPlaneWidget(x_coor=screen_x//2, y_coor=300, size=(400, 260), widget_color=(0, 100, 200))
    bkdisk_widget_one = BeltramiKleinDiskWidget(x_coor=screen_x//2, y_coor=300, size=(300, 300),color=(255,0,255))

    # setting up shapelist
    newshape = Hyperbolic_Shape()

    # making 2 default type 2 lines for testing purposes
    shapelist[selector].include(Hyperbolic_Line(complex(0.2, -0.4), complex(0.4, -0.2), tcolor=[0, 0, 255], tpath=False))
    shapelist[selector].include(Hyperbolic_Line(complex(0.5, -0.5), complex(0.6, -0.2), tpath=True))
    shapelist[selector].include(Hyperbolic_Line(complex(0, -0.5), complex(0, -0.2), tcolor=[0, 0, 255], tpath=False))

    #generate_n_random_lines(20)

    redraw_all()


def mousePressed():
    """
    Runs on any mouse click, position of the click is stored in "mouseX" and "mouseY"
    """
    global stored
    global clicks
    global selector
    global shapelist
    global mode
    global buttonlist
    if abs(interpret_point([mouseX, mouseY])) <= 1:
        if not isinstance(selector, int):
            draw_in = len(shapelist)
        else:
            draw_in = selector
        if clicks == 1:
            if mode <= 1:
                while len(shapelist) <= draw_in:
                    newshape = Hyperbolic_Shape()
                if mode == 0:  # draw line
                    shapelist[draw_in].include(
                        Hyperbolic_Line(stored, interpret_point([mouseX, mouseY]), [0, 0, 255], False))
                else:  # draw path
                    shapelist[draw_in].include(Hyperbolic_Line(stored, interpret_point([mouseX, mouseY])))
            elif mode == 3 and len(shapelist) >= selector:
                newshape = Hyperbolic_Shape()
                active_translation = P_disc.intuitive_translation_constructor(stored, interpret_point([mouseX, mouseY]))
                for l in shapelist[selector].lines:
                    newshape.include(l.intuitive_disk_translate(active_translation))
                shapelist.append(newshape)
            elif mode == 9:  # draw equal shape
                P_disc.equilateral_Shape(selector, abs(interpret_point([mouseX, mouseY]) - stored))
            clicks = 0
            redraw_all()
        else:
            stored = interpret_point([mouseX, mouseY])
            clicks += 1
    else:  # Button is (probably) meant to be pressed , check which one
        for button in buttonlist:
            if button.check([mouseX, mouseY]) == True:
                break
        for plane in plane_list:
            for button in plane.button_list:
                if button.check([mouseX, mouseY]) == True:
                    break
        for bkdisk_widget in bkdisk_list:
            for button in bkdisk_widget.button_list:
                if button.check([mouseX, mouseY]) == True:
                    break


def mouseReleased():
    for plane_widget in plane_list:
        delta = plane_widget.button_drag_drop.drop()
        if delta != (0, 0):
            plane_widget.drag(delta)
    for bkdisk_widget in bkdisk_list:
        delta = bkdisk_widget.button_drag_drop.drop()
        if delta != (0, 0):
            bkdisk_widget.drag(delta)


def keyPressed():
    """
    Runs on any keyboard press , pressed key is automatically stored in "key"
    """
    global input
    global typing
    global selector
    if key == ENTER:
        input = typing
        interpret_input(input)
        typing = ""
    else:
        typing = typing + str(key)
        clear_textbox()
        fill(0, 0, 0)
        text(str(typing), 0.6 * width / 20, 18.8 * height / 20, 3.5 * width / 20, height / 20)


def draw():
    """
    loop function, calls "pressed functions" on its own
    """
    return
