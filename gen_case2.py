
from steady_fit import FallingLimb_gen, RisingLimb_gen
from steady_fit import FallingLimb_fit, RisingLimb_fit, fit_steady


def run():
    """
    ks1 = [1.1, 1.5, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0, 100.0,
           120.0]
    ks2 = [1.1, 5.0, 10.0, 25.0, 50.0, 75.0, 90.0, 100.0, 105.0, 110.0, 115.0,
           120.0]
    # Rising Limb s-q generation and fit
    Rcase2 = RisingLimb_gen(10.0, 2.0, 0.1, 0.01, ks1, rec_name='ks2',
                            folder='Steady_fit', gen=True)
    # Falling Limb s-q generation and fit
    Fcase2 = FallingLimb_gen(10.0, 2.0, 0.1, 0.01, ks2, rec_name='ks2',
                             folder='Steady_fit', gen=True)
    Rcase2.run()
    Fcase2.run()
    """

    Rcase2_fit = RisingLimb_fit(rec_name='ks2', folder='Steady_fit',
                                save_name='ks2')
    Rcase2_fit.run()

    Fcase2_fit = FallingLimb_fit(rec_name='ks2', folder='Steady_fit',
                                 save_name='ks2')
    Fcase2_fit.run()

    ft = fit_steady('ks2', folder='Steady_fit')
    ft.param(10.0, 0.1, 2.0, 0.01)
    ft.steady_SQ()
