from statistics import median

INCIDENCE_ANGLES_ENVI = {
    'IMP': {
        'IS1': {
            'near' : 14.21,
            'far' : 22.37
        },
        'IS2': {
            'near' : 18.77,
            'far' : 26.18
        },
        'IS3': {
            'near' : 25.83,
            'far' : 31.2
        },
        'IS4': {
            'near' : 30.86,
            'far' : 36.14
        },
        'IS5': {
            'near' : 35.7,
            'far' : 39.22
        },
        'IS6': {
            'near' : 38.73,
            'far' : 42.63
        },
        'IS7': {
            'near' : 42.41,
            'far' : 45.19
        },
    },    
    'WSM' : {
        '-': {
            'near': 16.28,
            'far': 42.43
            }
    },
    'APP': {
        'IS1': {
            'near' : 14.23,
            'far' : 22.37
        },
        'IS2': {
            'near' : 20.77,
            'far' : 26.18
        },
        'IS3': {
            'near' : 27.43,
            'far' : 31.2
        },
        'IS4': {
            'near' : 31.61,
            'far' : 36.14
        },
        'IS5': {
            'near' : 36.82,
            'far' : 39.22
        },
        'IS6': {
            'near' : 39.74,
            'far' : 42.63
        },
        'IS7': {
            'near' : 43.25,
            'far' : 45.19
        },
    }
}

def get_incidence_angles(sensor, mode, swath_id):
    if sensor in ['ERS1', 'ERS2']:
        return 20.1, 25.9, 23
    if mode == 'WSM':
        swath_id = '-'
    near = INCIDENCE_ANGLES_ENVI[mode][swath_id]['near'] 
    far = INCIDENCE_ANGLES_ENVI[mode][swath_id]['far']
    mid = median([near, far])     
    return near, far, mid
    