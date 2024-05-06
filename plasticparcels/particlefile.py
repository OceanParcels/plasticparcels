from parcels import ParticleFile as pcls_particlefile


def ParticleFile(name, particleset, outputdt, settings, **kwargs):
    """Wrapper method to initialise a :class:`parcels.particlefile.ParticleFile` object from the ParticleSet."""
    pfile = pcls_particlefile(name=name, particleset=particleset, outputdt=outputdt, **kwargs)
    pfile.add_metadata('plasticparcels_settings', str(settings))
    return pfile
