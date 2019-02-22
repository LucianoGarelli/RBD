# -*- coding: utf-8 -*-
"""
@author: lgenzelis
"""

import numpy as np

class Quaternion:

        def __init__(self, d, v):
            """
            Crea el cuaternión a partir de la parte escalar y la parte vectorial
            """
            self.d = d
            self.v = np.asarray(v)

        def __str__(self):
            return '<{0}, ({1}, {2}, {3})>'.format(self.d, self.v[0], self.v[1], self.v[2])

        __repr__ = __str__

        @classmethod
        def fromAngleAndDirection(cls, alfa, u):
            """
            Crea el cuaternión a partir de un ángulo y un vector dirección unitario
            """
            return Quaternion(np.cos(alfa/2.), np.asarray(u) * np.sin(alfa/2.))

        @classmethod
        def fromRotationVector(cls, v):
            """
            Crea el cuaternión a partir de un vector de rotación
            """
            alfa = np.linalg.norm(v)
            if alfa == 0:
                u = np.zeros(3)
            else:
                u = v/alfa
            return cls.fromAngleAndDirection(alfa, u)

        @classmethod
        def fromEulerAngles(cls, roll, pitch, yaw):
            """
            Crea el cuaternión a partir de los ángulos de Euler (o más bien de Tait-Bryan). Se supone que la rotación correspondiente
            a los ángulos de Euler se realiza en el orden zyx. De esta forma, las evoluciones de mis modelos usando cuaterniones y
            ángulos de euler son equivalentes entre sí si los cuaterniones se obtienen a partir de los ángulos de Euler usando este
            constructor.
            """
            q_z_conj = Quaternion.fromAngleAndDirection(-yaw, [0, 0, -1])
            q_y_conj = Quaternion.fromAngleAndDirection(-pitch, [0, -1, 0])
            q_x_conj = Quaternion.fromAngleAndDirection(-roll, [-1, 0, 0])
            return q_z_conj.mult(q_y_conj).mult(q_x_conj)

        def mult(self, q, normalize=False):
            """
            Devuelve el producto del cuaternión por el cuaternión q
            :param q: cuaternión por el cual se postmultiplica
            :param normalize: si se normaliza el resultado. Si los cuaterniones son unitarios no haría falta (en teoría), pero por cuestiones numéricas podría ser aconsejable usarlo.
            :return: el producto
            """
            qr = Quaternion(self.d*q.d - np.dot(self.v, q.v), self.d*q.v + q.d*self.v + np.cross(self.v, q.v))
            if normalize:
                qr.normalize()
            return qr

        def modulus(self):
            return np.sqrt(self.d*self.d + np.dot(self.v, self.v))

        def normalize(self):
            modul = self.modulus()
            assert modul > 1e-8
            self.d /= modul
            self.v /= modul

        def get_conjugate(self):
            return Quaternion(self.d, -self.v)

        def calc_rot_matrix(self, unit_quat=True):
            """
            Calcula la matriz de rotación M tal que, dado un vector v y mi cuaternión q, es lo mismo hacer Mv que q*v*q_inv
            :param unit_quat: dice si el cuaternión es unitario o no
            """
            assert unit_quat
            # e1q = Quaternion(0., [1.,0,0])
            # e2q = Quaternion(0., [0,1.,0])
            # e3q = Quaternion(0., [0,0,1.])
            # c1 = self.mult(e1q).mult(self.get_conjugate())  # columna 1 de la matriz
            # c2 = self.mult(e2q).mult(self.get_conjugate())  # columna 2 de la matriz
            # c3 = self.mult(e3q).mult(self.get_conjugate())  # columna 3 de la matriz
            # Q = np.array([c1.v, c2.v, c3.v]).T

            #  de "The quaternions with an application to rigid body dynamics"
            q0 = self.d
            q1, q2, q3 = self.v
            q0_2 = q0**2
            q1_2 = q1**2
            q2_2 = q2**2
            q3_2 = q3**2
            q0q1 = q0*q1
            q0q2 = q0*q2
            q0q3 = q0*q3
            q1q2 = q1*q2
            q1q3 = q1*q3
            q2q3 = q2*q3
            # la obtención de Q me cuesta 16 multiplicaciones
            Q = np.array([[q0_2 + q1_2 - q2_2 - q3_2,   2.*(q1q2 - q0q3),            2.*(q1q3 + q0q2)],
                         [2.*(q1q2 + q0q3),             q0_2 - q1_2 + q2_2 - q3_2,   2.*(q2q3 - q0q1)],
                         [2.*(q1q3 - q0q2),             2.*(q2q3 + q0q1),            q0_2 - q1_2 - q2_2 + q3_2]])

            return Q

        def rotate_vector(self, vec, unit_quat=True):
            """
            Calcula la rotación del vector v mediante el cuaternión q a través de la fórmula q*v*q_inv
            :param unit_quat: dice si el cuaternión es unitario o no
            """
            #print  ('Warning: si se quiere hacer mas de una rotacion con el mismo cuaternion, es mucho mas rapido calcular la ' \
             #      'matriz de rotacion (calc_rot_matrix) y despues hacer el producto de la matriz por los vectores.')
            assert unit_quat
            q_vec = Quaternion(0, vec)
            return self.mult(q_vec).mult(self.get_conjugate()).v

        def get_angle_and_axys(self, unit_quat=True):
            assert unit_quat
            angle = np.arccos(self.d) * 2.  # esto esta entre 0 y 2pi
            if angle > np.pi:
                angle -= 2*np.pi  # para que el ángulo esté entre -pi y pi
            axys = self.v/np.linalg.norm(self.v)
            return angle, axys

        def get_euler_anles(self, order='zyx'):
            """
            Devuelve los ángulos de Euler, tal que, si se realiza una rotacion con dichos ángulos dé el mismo resultado que si la rotación
            se realiza con el cuaternión. Para que el resultado sea único, acá se supone que el roll y el yaw están en
            [-pi, pi], y el pitch en [-pi/2, pi/2]. Si se obtuviese el cuaternión usando ángulos de Euler fuera de estos rangos, y luego
            se llamase a esta función, los ángulos resultado no coincidirían con los originales.
            :param order: el orden en que se realizaría la rotación con los ángulos de Euler (para dar el mismo resultado
            que realizar la rotación con el cuaternión). Puede ser 'zyx' o 'xyz'.
            :return: Los ángulos de euler (o más bien de Tait-Bryan): roll (x), pitch (y), yaw (z)
            """
            if order == 'zyx':
                q0 = self.d
                q1, q2, q3 = self.v
                roll = np.arctan2(q2*q3+q0*q1, .5 - (q1**2+q2**2))
                pitch = np.arcsin(-2.*(q1*q3-q0*q2))
                yaw = np.arctan2(q1*q2+q0*q3, .5 - (q2**2+q3**2))
            elif order == 'xyz':
                q0 = self.d
                q1, q2, q3 = -self.v
                roll = -np.arctan2(q2*q3+q0*q1, .5 - (q1**2+q2**2))
                pitch = -np.arcsin(-2.*(q1*q3-q0*q2))
                yaw = -np.arctan2(q1*q2+q0*q3, .5 - (q2**2+q3**2))
            else:
                raise Exception('El orden "{0}" es desconocido'.format(order))
            return roll, pitch, yaw

        def __mul__(self, escalar):
            assert isinstance(escalar, (int,float))
            return Quaternion(self.d*escalar, self.v*escalar)

        __rmul__ = __mul__

        def mult_cuat_times_vec(self, v):
            """
            Multiplica al cuaternion por el vector v (en ese orden). Es equivalente a self.mult(Quaternion(0, v)), pero es mas rapido.
            :param v: vector en R3
            :return: el cuaternion resultado de multiplicar este cuaternion por el cuaternion <0, v>
            """
            q0 = self.d
            q1, q2, q3 = self.v
            q0_r = -v[0]*q1 - v[1]*q2 - v[2]*q3
            q1_r = v[0]*q0 + v[2]*q2 - v[1]*q3
            q2_r = v[1]*q0 - v[2]*q1 + v[0]*q3
            q3_r = v[2]*q0 + v[1]*q1 - v[0]*q2
            return Quaternion(q0_r, [q1_r, q2_r, q3_r])

