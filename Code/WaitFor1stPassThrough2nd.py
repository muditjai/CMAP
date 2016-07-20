import sys

if __name__ == '__main__':
    print(sys.argv)
    with open(sys.argv[2]) as fp:
        with open(sys.argv[3], mode='w') as write_fp:
            write_fp.write(fp.read(-1))
