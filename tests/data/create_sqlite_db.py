"""
This script is intended to create SQLite databases for use in legacy import tests.
"""
import os
from pyroSAR import Archive
from spatialist.ancillary import finder

import logging

logging.basicConfig(level=logging.DEBUG)

target = os.path.dirname(os.path.realpath(__file__))

scenes = finder(
    target=target,
    matchlist=['S1A*.zip', 'ASA*', 'SAR*', '00000*.zip'],
    recursive=False
)

out = os.path.join(target, 'archive_outdated_tables-unmanaged.db')

with Archive(out) as db:
    db.insert(scene_in=scenes)
    assert db.size == (8, 1)
    
    # Make the fixture portable by storing scene paths relative to the
    # directory containing the database/test data.
    with db.engine.begin() as conn:
        for table in ['data', 'duplicates']:
            rows = conn.exec_driver_sql(
                f'SELECT scene FROM {table}'
            ).fetchall()
            
            for row in rows:
                scene = row[0]
                relative = os.path.relpath(scene, start=target)
                
                conn.exec_driver_sql(
                    f'UPDATE {table} SET scene = ? WHERE scene = ?',
                    (relative, scene),
                )
