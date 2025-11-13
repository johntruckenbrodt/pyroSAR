from typing import Literal

RESOLUTION_NESZ = {
    'ERS1': {
        'IMP': {
            'IS2': {
                'res_rg': 25.04,
                'res_az': 21.51,
                'nesz_nr': 26.8,
                'nesz_fr': 26
            },
            'std_dev': 20
        },
        'IMS': {
            'IS2': {
                'res_rg': 5.32,
                'res_az': 9.66,
                'nesz_nr': 26.8,
                'nesz_fr': 26
            },
            'std_dev': 20
        }
    },
    'ERS2': {
        'IMP': {
            'IS2': {
                'res_rg': 21.63,
                'res_az': 25.19,
                'nesz_nr': 23.1,
                'nesz_fr': 21.5
            },
            'std_dev': 20
        },
        'IMS': {
            'IS2': {
                'res_rg': 5.33,
                'res_az': 9.83,
                'nesz_nr': 23.1,
                'nesz_fr': 21.5
            },
            'std_dev': 20
        }
    },
    'ASAR': {
        'IMP': {
            'IS1': {
                'res_rg': 30.86,
                'res_az': 22.14,
                'nesz_nr': 25.1,
                'nesz_fr': 19.2
            },
            'IS2': {
                'res_rg': 24.90,
                'res_az': 22.14,
                'nesz_nr': 21.8,
                'nesz_fr': 20.5
            },
            'IS3': {
                'res_rg': 24.84,
                'res_az': 22.14,
                'nesz_nr': 22.6,
                'nesz_fr': 20.5
            },
            'IS4': {
                'res_rg': 25.56,
                'res_az': 22.14,
                'nesz_nr': 22.3,
                'nesz_fr': 19.1
            },
            'IS5': {
                'res_rg': 25.73,
                'res_az': 22.14,
                'nesz_nr': 21.4,
                'nesz_fr': 19
            },
            'IS6': {
                'res_rg': 26.15,
                'res_az': 22.14,
                'nesz_nr': 24,
                'nesz_fr': 21.2
            },
            'IS7': {
                'res_rg': 26.59,
                'res_az': 22.14,
                'nesz_nr': 23,
                'nesz_fr': 20.4
            },
            'std_dev': 5
        },
        'IMS': {
            'IS1': {
                'res_rg': 5.77,
                'res_az': 8.43,
                'nesz_nr': 25.1,
                'nesz_fr': 19.2
            },
            'IS2': {
                'res_rg': 5.77,
                'res_az': 8.43,
                'nesz_nr': 21.8,
                'nesz_fr': 20.5
            },
            'IS3': {
                'res_rg': 5.77,
                'res_az': 8.43,
                'nesz_nr': 22.6,
                'nesz_fr': 20.5
            },
            'IS4': {
                'res_rg': 5.77,
                'res_az': 8.43,
                'nesz_nr': 22.3,
                'nesz_fr': 19.1
            },
            'IS5': {
                'res_rg': 5.77,
                'res_az': 8.43,
                'nesz_nr': 21.4,
                'nesz_fr': 19
            },
            'IS6': {
                'res_rg': 5.77,
                'res_az': 8.43,
                'nesz_nr': 24,
                'nesz_fr': 21.2
            },
            'IS7': {
                'res_rg': 5.77,
                'res_az': 8.43,
                'nesz_nr': 23,
                'nesz_fr': 20.4
            },
            'std_dev': 5
        },
        'APP': {
            'IS1': {
                'res_rg': 31.22,
                'res_rg_new': 31.22,
                'res_az': 27.45,
                'nesz_nr': 25.1,
                'nesz_fr': 19.2
            },
            'IS2': {
                'res_rg': 25.23,
                'res_rg_new': 24.10,
                'res_az': 27.45,
                'nesz_nr': 21.8,
                'nesz_fr': 20.5
            },
            'IS3': {
                'res_rg': 24.74,
                'res_rg_new': 24.30,
                'res_az': 27.45,
                'nesz_nr': 22.6,
                'nesz_fr': 20.5
            },
            'IS4': {
                'res_rg': 25.46,
                'res_rg_new': 25.30,
                'res_az': 27.45,
                'nesz_nr': 22.3,
                'nesz_fr': 19.1
            },
            'IS5': {
                'res_rg': 25.70,
                'res_rg_new': 25.35,
                'res_az': 27.45,
                'nesz_nr': 21.4,
                'nesz_fr': 19
            },
            'IS6': {
                'res_rg': 26.07,
                'res_rg_new': 25.90,
                'res_az': 27.45,
                'nesz_nr': 24,
                'nesz_fr': 21.2
            },
            'IS7': {
                'res_rg': 26.53,
                'res_rg_new': 26.32,
                'res_az': 27.45,
                'nesz_nr': 23,
                'nesz_fr': 20.4
            },
            'std_dev': 10
        },
        'APS': {
            'IS1': {
                'res_rg': 4.3,
                'res_az': 8.39,
                'nesz_nr': 25.1,
                'nesz_fr': 19.2
            },
            'IS2': {
                'res_rg': 4.3,
                'res_az': 8.39,
                'nesz_nr': 21.8,
                'nesz_fr': 20.5
            },
            'IS3': {
                'res_rg': 4.3,
                'res_az': 8.39,
                'nesz_nr': 22.6,
                'nesz_fr': 20.5
            },
            'IS4': {
                'res_rg': 4.3,
                'res_az': 8.39,
                'nesz_nr': 22.3,
                'nesz_fr': 19.1
            },
            'IS5': {
                'res_rg': 4.3,
                'res_az': 8.39,
                'nesz_nr': 21.4,
                'nesz_fr': 19
            },
            'IS6': {
                'res_rg': 4.3,
                'res_az': 8.39,
                'nesz_nr': 24,
                'nesz_fr': 21.2
            },
            'IS7': {
                'res_rg': 4.3,
                'res_az': 8.39,
                'nesz_nr': 23,
                'nesz_fr': 20.4
            },
            'std_dev': 10
        },
        'WSM': {
            'WS': {
                'res_rg': 150,
                'res_az': 150,
                'nesz_nr': 19.5,
                'nesz_fr': 23.5
            },
            'std_dev': 20
        },
        'WSS': {
            'WS': {
                'res_rg': None,
                'res_az': None,
                'nesz_nr': None,
                'nesz_fr': None
            },
            'std_dev': None
        }
    }
}


def get_resolution_nesz(
        sensor: Literal['ERS1', 'ERS2', 'ASAR'],
        mode: Literal['APP', 'APS', 'IMP', 'IMS', 'WSM', 'WSS'],
        swath_id: Literal['IS1', 'IS2', 'IS3', 'IS4', 'IS5', 'IS6', 'IS7', 'WS'],
        date: str
) -> tuple[int | float | None, int | float | None, int | float | None, int | float | None]:
    """
    Get acquisition characteristics not contained in the product metadata:

    - range resolution
    - azimuth resolution
    - near range noise equivalent sigma zero (NESZ)
    - far range NESZ
    
    Parameters
    ----------
    sensor:
        the satellite sensor
    mode:
        the sensor acquisition mode
    swath_id:
        the sensor swath ID
    date:
        the acquisition date formatted as YYYYmmdd/YYYYmmddTHHMMSS

    Returns
    -------
        the attributes listed above
    """
    suffix = '_new' if mode == 'APP' and date > '20090528' else ''
    data = RESOLUTION_NESZ[sensor][mode][swath_id]
    return (data[f'res_rg{suffix}'], data['res_az'],
            data['nesz_nr'], data['nesz_fr'])
