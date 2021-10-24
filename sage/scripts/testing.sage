
P.<s_13_x, s_13_y, t_13_x, t_13_y, d_13_x, d_13_y, u_13_1, u_13_2, t_14_x, t_14_y> = PolynomialRing(RR, 10)

s_13_x, s_13_y, t_13_x, t_13_y, d_13_x, d_13_y, u_13_1, u_13_2, t_14_x, t_14_y = var('s_13_x s_13_y t_13_x t_13_y d_13_x d_13_y u_13_1 u_13_2 t_14_x t_14_y')

poly_0 = t_13_x + -1 * s_13_x
poly_1 = t_13_y + -1 * s_13_y
poly_2 = -1 * t_14_x + d_13_x*u_13_2 + s_13_x + -1 * s_13_x*u_13_2
poly_3 = -1 * t_14_y + d_13_y*u_13_2 + s_13_y + -1 * s_13_y*u_13_2
poly_4 = -5.242502821306275 + -0.873911842479982 * u_13_1 + -1 * d_13_x + 0.7142857142857143 * s_13_y + s_13_x + 0.2857142857142857 * s_13_x*u_13_1
poly_5 = -3.276185529910945 + -2.554455308417708 * u_13_1 + -1 * d_13_y + s_13_y + 0.2857142857142857 * s_13_y*u_13_1 + -0.7142857142857143 * s_13_x
poly_6 = -63.45673860828117 * u_13_1 + 89.28980693134451 * u_13_1^2 + 7.64535119055526 * s_13_y*u_13_1 + -17.88118715892396 * s_13_y*u_13_1^2 + s_13_y^2*u_13_1^2 + -1.601089629633194 * s_13_x*u_13_1 + -6.117382897359875 * s_13_x*u_13_1^2 + s_13_x^2*u_13_1^2
poly_7 = 68.9530141674729 + -14.67900789965757 * s_13_y + s_13_y^2 + 9.173319483750646 * s_13_x + s_13_x^2
poly_8 = -10.03592973490771 * t_14_y + t_14_y^2 + -13.69309700325307 * t_14_x + t_14_x^2 + 10.03592973490771 * s_13_y + -1 * s_13_y*t_14_y + 13.69309700325307 * s_13_x + -1 * s_13_x*t_14_x
poly_9 = 285.7085753283435 + -20.07185946981542 * t_14_y + t_14_y^2 + -27.38619400650613 * t_14_x + t_14_x^2



P.<s_13_x, s_13_y, t_13_x, t_13_y, d_13_x, d_13_y, u_13_1, u_13_2, t_14_x, t_14_y> = PolynomialRing(QQ, 10)

sols_ideal = P.ideal([poly_0, poly_1, poly_2, poly_3, poly_4, poly_5, poly_6, poly_7, poly_8, poly_9])
groebner = sols_ideal.groebner_basis()

sols_ideal.dimension()
sols_ideal.variety()


solve(
	[poly_0, poly_1, poly_2, poly_3, poly_4, poly_5, poly_6, poly_7, poly_8, poly_9],
	s_13_x, s_13_y, t_13_x, t_13_y, d_13_x, d_13_y, u_13_1, u_13_2, t_14_x, t_14_y
)



P.<s_13_x, s_13_y, t_13_x, t_13_y, d_13_x, d_13_y, u_13_1, u_13_2, t_14_x, t_14_y> = PolynomialRing(RR, 10)

ax = -20.0; ay = 20.0;
px = 0.0; py = -20.0;
r1 = 2.26 / 2;
r2 = 2.26 / 2;
alpha = 2.0 / 7;

tan_dir_x = y;
tan_dir_y = -x;

tx = x + tan_dir_x;
ty = y + tan_dir_y;

aim_dir_x = x - ax;
aim_dir_y = y - ay;

aim_x = x + s1 * (aim_dir_x - x);
aim_y = y + s1 * (aim_dir_y - y);

dest_x = x + s2 * (dx - x);
dest_y = y + s2 * (dy - y);

c_glance_1 = (tx - x) * alpha + (aim_x - y) * (1 - alpha) - (dx - x);
c_glance_2 = (ty - y) * alpha + (aim_y - y) * (1 - alpha) - (dy - y);
c_orth_req_1 = (aim_x - x) * (aim_x - tx) + (aim_y - y) * (aim_y - ty);
c_orth_req_2 = (dest_x - x) * (dest_x - px) + (dest_y - y) * (dest_y - py);
c_radius_1 = x^(2) +y^(2) - (r1 + r2)^2;
c_radius_2 = (dest_x - px)^(2) + (dest_y - py)^2 - r1^2;


c_glance_1 - glance_1
c_glance_2 - glance_2
c_orth_req_1 - orth_req_1
c_orth_req_2 - orth_req_2
c_radius_1 - radius_1
c_radius_2 - radius_2

poly_0 = -1 * s_13_x + t_13_x
poly_1 = -1 * s_13_y + t_13_y
poly_2 = -1 * s_13_x * u_13_2 + s_13_x + d_13_x * u_13_2 + -1 * t_14_x
poly_3 = -1 * s_13_y * u_13_2 + s_13_y + d_13_y * u_13_2 + -1 * t_14_y
poly_4 = 0.2857142857142857 * s_13_x * u_13_1 + s_13_x + 0.7142857142857143 * s_13_y + -1 * d_13_x + 0.4314528043649248 * u_13_1 + -5.278825623728934
poly_5 = -0.7142857142857143 * s_13_x + 0.2857142857142857 * s_13_y * u_13_1 + s_13_y + -1 * d_13_y + 1.337511970443132 * u_13_1 + 4.639472547824646
poly_6 = s_13_x^2 * u_13_1^2 + 3.020169630554474 * s_13_x * u_13_1^2 + 12.07164776977147 * s_13_x * u_13_1 + s_13_y^2 * u_13_1^2 + 9.362583793101928 * s_13_y * u_13_1^2 + -8.005346382231741 * s_13_y * u_13_1 + 24.19484997004461 * u_13_1^2 + -19.2461511557178 * u_13_1
poly_7 = s_13_x^2 + -12.99052313390901 * s_13_x + s_13_y^2 + -14.78071174644101 * s_13_y + 96.51877361842173
poly_8 = -1 * s_13_x * t_14_x + 6.272184894087336 * s_13_x + -1 * s_13_y * t_14_y + 7.115132186973382 * s_13_y + t_14_x^2 + -6.272184894087336 * t_14_x + t_14_y^2 + -7.115132186973382 * t_14_y
poly_9 = t_14_x^2 + -12.54436978817467 * t_14_x + t_14_y^2 + -14.23026437394676 * t_14_y + 89.82527074891658



sols_ideal = P.ideal([poly_0, poly_1, poly_2, poly_3, poly_4, poly_5, poly_6, poly_7, poly_8, poly_9])
sols_ideal.dimension()
sols_ideal.variety()
