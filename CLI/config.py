from importlib.metadata import distributions

for dist in distributions():
    print(dist.metadata['Name'], dist.metadata['Version'])
