import VOEvent-ReaderTst as vread
import os
import voeventparse
import sys

def main():
    stdin = sys.stdin.read()
    v = voeventparse.loads(stdin)
    vread.handle_voevent(v)
    return 0
    
    
if __name__ == '__main__':
    sys.exit(main())