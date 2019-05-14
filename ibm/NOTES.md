# Notes sur STAMM en vue de la refactorisation

* Récupérer la date de début de la simulation à partir des fichiers ou des paramètres
* Clarifier la façon dont stamm boucle après le dernier fichier de forçage disponible:

Vérifier comment cette condition (dans `IOlib.py`) est géré dans le module d'avection:
```{python}
if param['nsteps_simu']+np.max(t_init)>param['nsteps_max']:
        print 'WARNING - There are not enough input files : Loop on time when last one is reached'
        print 'WARNING - max(t_init) = %d ' % np.max(t_init)
        print 'WARNING - nsteps_max  = %d ' % (param['nsteps_max'])
        print 'WARNING - nsteps_simu = %d, should be less than %d ' % (param['nsteps_simu'],param['nsteps_max']-np.max(t_init))
```

Dans la namelist, je pense qu'il vaut mieux spécifier nsteps_max=nombre de fichiers de forçage disponibles et adapter nsteps_simu en foncton de max(t_init)

