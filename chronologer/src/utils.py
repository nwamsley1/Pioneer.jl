'''
Basic utility functions
'''

import time, datetime


def timer( start, ):
    return str( datetime.timedelta( seconds=round( time.time() - start ) ) )

def is_float( element ):
    try:
        float( element )
        return True
    except ValueError:
        return False