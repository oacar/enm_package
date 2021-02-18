import pytest
from src.enm import Enm
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

@pytest.fixture
def enm():
    """ Returns empty Enm object with name Empty"""
    enm =  Enm('Empty')

    L = np.array([[ 5,  0,  0,  0,  0, -1, -1, -1, -1, -1],
        [ 0,  5,  0,  0,  0, -1, -1, -1, -1, -1],
        [ 0,  0,  3,  0,  0, -1, -1,  0,  0, -1],
        [ 0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
        [ 0,  0,  0,  0,  4, -1, -1, -1,  0, -1],
        [-1, -1, -1, -1, -1,  8, -1, -1, -1,  0],
        [-1, -1, -1,  0, -1, -1,  7, -1, -1,  0],
        [-1, -1,  0,  0, -1, -1, -1,  6, -1,  0],
        [-1, -1,  0,  0,  0, -1, -1, -1,  6, -1],
        [-1, -1, -1,  0, -1,  0,  0,  0, -1,  5]])
    enm.L=L
#L = D-A
#A = D-L
    adj = np.diag(np.diag(L)) - L
    enm.G = nx.from_numpy_array(adj)
    enm.giant_component()
    return enm


# def test_read_network(enm):
#     fname, path = enm.read_network('this/random/file/path.csv')

#     assert fname == 'this/random/file/path'
#     assert path =='.csv'


from unittest.mock import Mock
#from .. import my_module
def test_gnm(enm):
    enm.get_gnm()
    eigvals = enm.gnm.getEigvals()
    #eigvecs = enm.gnm.getEigvecs()
    np.testing.assert_array_almost_equal(eigvals,np.array([0.98062692, 2.93886417, 4.        , 5.        , 5.01742587,
       7.        , 7.38411109, 8.37647156, 9.30250039]))

def test_prs(enm):
    enm.get_gnm()
    enm.get_prs()
    prs_mat_got = enm.prs_mat
    prs_mat_expected = np.array([[0.00000000e+00, 1.00143605e-02, 3.90208376e-02, 3.75411396e-01,
        1.94212248e-02, 3.92782451e-03, 3.69922763e-05, 1.23696857e-03,
        3.25827872e-03, 1.96429686e-04],
       [1.00143605e-02, 0.00000000e+00, 3.90208376e-02, 3.75411396e-01,
        1.94212248e-02, 3.92782451e-03, 3.69922763e-05, 1.23696857e-03,
        3.25827872e-03, 1.96429686e-04],
       [1.43630592e-02, 1.43630592e-02, 0.00000000e+00, 1.40579881e-01,
        1.43630592e-02, 1.70009031e-03, 2.31658508e-07, 2.02513037e-02,
        1.43630592e-02, 1.39738986e-03],
       [1.52005688e-02, 1.52005688e-02, 1.54641195e-02, 0.00000000e+00,
        1.52005688e-02, 1.50887580e-05, 1.45516022e-02, 1.44235077e-02,
        1.52005688e-02, 2.12108155e-02],
       [1.27183177e-02, 1.27183177e-02, 2.55534558e-02, 2.45844505e-01,
        0.00000000e+00, 2.57220235e-03, 2.42250181e-05, 8.10049801e-04,
        1.27183177e-02, 1.28635304e-04],
       [1.21173773e-02, 1.21173773e-02, 1.42487966e-02, 1.14962887e-03,
        1.21173773e-02, 0.00000000e+00, 7.54388663e-03, 6.75862691e-03,
        1.21173773e-02, 9.31199383e-02],
       [7.40385244e-05, 7.40385244e-05, 1.25963463e-06, 7.19291822e-01,
        7.40385244e-05, 4.89424038e-03, 0.00000000e+00, 9.64180330e-04,
        7.40385244e-05, 2.20588616e-02],
       [1.73502383e-03, 1.73502383e-03, 7.71701603e-02, 4.99649297e-01,
        1.73502383e-03, 3.07290190e-03, 6.75706891e-04, 0.00000000e+00,
        1.73502383e-03, 2.32626981e-02],
       [4.58659989e-03, 4.58659989e-03, 5.49286862e-02, 5.28457511e-01,
        2.73387868e-02, 5.52910325e-03, 5.20731296e-05, 1.74125064e-03,
        0.00000000e+00, 2.76509303e-04],
       [1.91193930e-04, 1.91193930e-04, 3.69516813e-03, 5.09884672e-01,
        1.91193930e-04, 2.93800886e-02, 1.07276198e-02, 1.61428527e-02,
        1.91193930e-04, 0.00000000e+00]])
    np.testing.assert_array_almost_equal(prs_mat_expected,prs_mat_got,err_msg='PRS matrix is wrong')


def test_create_df(enm):

    enm.get_gnm()
    enm.get_prs()
    enm.create_df()
    got = enm.df.deg
    asked =  enm.degree   #assert got.shape==asked.shape
   # print(enm.df['smallest_eigenvec'])
    np.testing.assert_array_almost_equal(np.abs(got),np.abs(asked))


# def test_plot_network_spring(enm):
#     adj = np.array([[0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
#         [0, 0, 0, 0, 0, 1, 1, 0, 1, 0],
#         [0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
#         [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
#         [0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
#         [1, 1, 1, 1, 1, 0, 1, 0, 1, 1],
#         [1, 1, 1, 1, 0, 1, 0, 1, 1, 0],
#         [1, 0, 1, 1, 1, 0, 1, 0, 0, 1],
#         [1, 1, 0, 1, 0, 1, 1, 0, 0, 1],
#         [1, 0, 0, 1, 0, 1, 0, 1, 1, 0]])
#     enm.G = nx.from_numpy_array(adj)
#     enm.giant_component()
#     plot_network_spring

def test_get_sensor_effector(enm):
  enm.get_sensor_effector(use_threshold=True, quantile_threshold=0.99)
  sensors_df = np.array([[ 3.00000000e+00,  1.00000000e+00,  5.32351458e-02,
        -9.94778394e-02,  2.16280074e-02,  5.22439257e-02,
        -9.02056208e-16,  5.75928194e-16,  7.51687882e-03,
        -3.29597460e-15,  1.26467409e-01,  3.39568011e+00,
         0.00000000e+00,  0.00000000e+00,  7.35503582e-02,
         5.00000000e-01, -9.40221504e-01]])
  effectors_df = np.array([[ 6.00000000e+00,  7.00000000e+00,  7.91791550e-01,
         1.15440817e-01, -3.97381472e-01, -2.97773277e-01,
         3.05311332e-15, -2.49800181e-16, -1.89390545e-02,
         1.48769885e-14,  7.47506518e-01,  3.36493299e-02,
         7.73148148e-02,  5.71428571e-01,  4.02157939e-01,
         8.18181818e-01,  1.13136700e-01]])
  np.testing.assert_array_almost_equal(sensors_df,enm.sensors_df.values)
  np.testing.assert_array_almost_equal(effectors_df,enm.effectors_df.values)

def is_dataframe_big(dataframe):

    # Lets check if this dataframe has over a million rows
    if dataframe.shape[0] > 1000000:
        return True
    else:
        return False

def test_big():
    dataframe = Mock()
    dataframe.shape = (2000000, 1)

    assert is_dataframe_big(dataframe) is True


def test_not_big():
    dataframe = Mock()
    dataframe.shape = (100, 1)

    assert is_dataframe_big(dataframe) is False