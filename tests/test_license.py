import os
import re
from datetime import datetime

def test_license_year():
    dir_current = os.path.dirname(os.path.abspath(__file__))
    license = os.path.join(dir_current, '..', 'LICENSE.txt')
    assert os.path.isfile(license)
    with open(license, 'r') as f:
        content = f.read()
        start, end = re.search('([0-9]{4})-([0-9]{4})', content).groups()
    year_current = datetime.now().year
    assert int(start) == 2014
    assert int(end) == year_current
