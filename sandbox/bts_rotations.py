
def main():
    bts_keys = [256, 2, 768, 3, 6, 8192, 2048, 4096, 31744,4, 32760, 28672, 1024, 32736, 96, 32764, 3072, 7,224,32256, 32672, 8, 32640, 32744, 12288, 31232, 384, 1536,128, 16, 160, 30720, 24, 32, 24576, 6144, 32704,32752, 32000, 32512, 512, 64, 5, 192, 20480]
    bts_keys.sort()
    print(bts_keys)
    print('size:', len(bts_keys))
    step = 4
    print('step:', [bts_keys[i] for i in range(0,len(bts_keys),step)])
    max = 0
    deltas = []
    for i in range(len(bts_keys)-1):
        delta = bts_keys[i+1]- bts_keys[i]
        deltas.append(delta)
        if delta > max:
            max = delta
    print("max:", max)
    print("deltas", deltas)

    print("uniq.deltas", sorted(set(deltas)))

if __name__ == '__main__':
    main()
