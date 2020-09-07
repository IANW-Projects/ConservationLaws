def generate_settings(map, keys):
    # based on the matlab counterpart
    # Atomatically generates a formatted settings string that contains all
    # compiler optimizations and OpenCL defines of map with the given keys where 
    # keys is a list, e.g. settings_mesh = generate_settings(I_Mesh, ['DX', 'DY', 'DZ'])

    settings = '';
    for key in keys:
        if key in map:
            if not(type(map[key]) is float or type(map[key]) is int):
                if key==  'optimizations':
                    settings = settings + map[key]
                elif 'USE' in map[key]:
                    settings = settings + ' -D' +  map[key] + '=1'
                else:
                    settings = settings + ' -D' + key + '=' + str(map[key])
            else:
                settings = settings + ' -D' + key + "=" + str(map[key])
        else:
            print('Unknown identifier ' +  key)
    return settings


