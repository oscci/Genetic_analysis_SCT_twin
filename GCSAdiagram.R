library(DiagrammeR)
library(DiagrammeRsvg)
#x and y coords for SNP rectangles are hand-crafted by trial and error
#quartz()
grViz("
      digraph SEM {
      
      graph [layout = neato,
      overlap = true,
      outputorder = edgesfirst]
      
      node [shape = rectangle,
      fontname = Helvetica]
      
      a [pos = '0,7!', label = 'SNP_1']
      b [pos = '-1,7!', label = 'SNP_2']
      c [pos = '-1.9,7!', label = 'SNP_3']
      d [pos = '-2.8,7!', label = 'SNP_4']
      e [pos = '-3.5,6.3!', label = 'SNP_5']
      f [pos = '-3.9,5.5!', label = 'SNP_6']
      g [pos = '-4.3,4.7!', label = 'SNP_7']
      h [pos = '-4.6,3.8!', label = 'SNP_8']
      i [pos = '-5,3!', label = 'SNP_9']
      j [pos = '-5.1,2!', label = 'SNP_10']
      k [pos = '-5.2,1!', label = 'SNP_11']
      l [pos = '-5.3,0!', label = 'SNP_12']
      m [pos = '-5.2,-1!', label = 'SNP_13']
      n [pos = '-5.1,-2!', label = 'SNP_14']
      o [pos = '-5,-3!', label = 'SNP_15']
      p [pos = '-4.6,-3.8!', label = 'SNP_16']
      q [pos = '-4.3,-4.7!', label = 'SNP_17']
      r [pos = '-3.8,-5.6!', label = 'SNP_18']
      s [pos = '-3,-6.2!', label = 'SNP_19']
      t [pos = '-2,-6.7!', label = 'SNP_20']
      u [pos = '-1,-6.7!', label = 'SNP_21']
      v [pos = '0,-6.7!', label = 'SNP_22']
 
      w [pos = '-1,0!', label = 'Gene', shape = ellipse,fontsize=20]
      x [pos = '2,0!', label = 'Neuro factor', shape = ellipse,fontsize=20]

      y [pos = '4.5,2!', label = 'PhenoA',fontsize=18]
      z [pos = '4.5,0!', label = 'PhenoB',fontsize=18]
      aa [pos = '4.5,-2!', label = 'PhenoC',fontsize=18]

      
      w->x
      w->a 
      w->b 
      w->c 
       w->d
      w->e
      w->f
      w->g 
      w->h
      w->i 
      w->j 
      w->k 
      w->l 
      w->m
      w->n
      w->o
      w->p
      w->q 
      w->r 
      w->s
      w->t
      w->u
      w->v
      x->y 
      x->z 
      x->aa 
   
      
      }
      ")

#NB problems with exporting diagram from GraphViz; currently using screenshot

