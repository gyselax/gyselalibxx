# Polar Splines

This folder contains methods specific to the manipulation of polar splines. The classes in this folder are analogous to the spline methods in DDC. Additionally classes are provided to represent the spline itself. These classes are data storage classes only. They group two fields. The first describes the coefficients in front of the B-splines which traverse the O-point, the second describes the coefficients in front of the other B-splines. This separation is used as the coefficients in front of the B-splines which don't traverse the O-point are more useful in a 2D field.
