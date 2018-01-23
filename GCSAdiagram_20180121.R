library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)
#Updated 23 Jan 2018 for final N SNPs
#x and y coords for SNP rectangles are hand-crafted by trial and error
#quartz()
my_graph<-grViz("
      digraph SEM {
      
      graph [layout = neato,
      overlap = true,
      outputorder = edgesfirst]
      
      node [shape = rectangle,
      fontname = Helvetica]
      
      a [pos = '-0,5.6!', label = 'SNP_25']
      b [pos = '1,5.6!', label = 'SNP_26']
      c [pos = '1.9,5.6!', label = 'SNP_27']
      d [pos = '2.8,5.5!', label = 'SNP_28']
      e [pos = '3.5,4.9!', label = 'SNP_29']
      f [pos = '3.9,4.2!', label = 'SNP_30']
      g [pos = '4.3,3.5!', label = 'SNP_31']
      h [pos = '4.6,2.8!', label = 'SNP_32']
      i [pos = '5,2.1!', label = 'SNP_33']
      j [pos = '5.1,1.4!', label = 'SNP_34']
      k [pos = '5.2,.7!', label = 'SNP_35']
      l [pos = '5.3,0!', label = 'SNP_36']
      m [pos = '5.2,-.7!', label = 'SNP_37']
      n [pos = '5.1,-1.4!', label = 'SNP_38']
      o [pos = '5,-2.1!', label = 'SNP_39']
      p [pos = '4.6,-2.7!', label = 'SNP_40']
      q [pos = '4.3,-3.3!', label = 'SNP_41']
      r [pos = '3.8,-3.9!', label = 'SNP_42']
      s [pos = '3.4,-4.5!', label = 'SNP_43']
      t [pos = '2.5,-5!', label = 'SNP_44']
      u [pos = '1.5,-5!', label = 'SNP_45']
      v [pos = '.5,-5!', label = 'SNP_46']
      vv [pos = '-0.5,-5!', label = 'SNP_47']
 
      w [pos = '-1,0!', label = 'NRXN1', shape = ellipse,fontsize=20]
      aw [pos = '-6,0!', label = 'CNTNAP2', shape = ellipse,fontsize=20]
      x [pos = '-3.3,-.5!', label = 'Neurodev\n factor', shape = ellipse,fontsize=20]

      y [pos = '-4.8,-2!', label = 'Nonword\nrepetition',fontsize=18]
      z [pos = '-3.3,-2!', label = 'Language\nfactor',fontsize=18]
      az [pos = '-1.8,-2!', label = 'Global\nimpairment',fontsize=18]

      a1 [pos = '-4,5.6!', label = 'SNP_1']
      aa [pos = '-5,5.6!', label = 'SNP_2']
      ab [pos = '-6,5.6!', label = 'SNP_3']
      ac [pos = '-6.9,5.6!', label = 'SNP_4']
      ad [pos = '-7.9,5.5!', label = 'SNP_5']
      ae [pos = '-8.5,4.9!', label = 'SNP_6']
      af [pos = '-8.9,4.2!', label = 'SNP_7']
      ag [pos = '-9.3,3.5!', label = 'SNP_8']
      ah [pos = '-9.6,2.8!', label = 'SNP_9']
      ai [pos = '-10,2.1!', label = 'SNP_10']
      aj [pos = '-10.1,1.4!', label = 'SNP_11']
      ak [pos = '-10.2,.7!', label = 'SNP_12']
      al [pos = '-10.3,0!', label = 'SNP_13']
      am [pos = '-10.2,-.7!', label = 'SNP_14']
      an [pos = '-10.1,-1.4!', label = 'SNP_15']
      ao [pos = '-10,-2.1!', label = 'SNP_16']
      ap [pos = '-9.8,-2.7!', label = 'SNP_17']
      aq [pos = '-9.6,-3.3!', label = 'SNP_18']
      ar [pos = '-9.4,-3.9!', label = 'SNP_19']
      a3 [pos = '-9,-4.5!', label = 'SNP_20']
      as [pos = '-8,-5!', label = 'SNP_21']
      at [pos = '-7,-5!', label = 'SNP_22']
      au [pos = '-6,-5!', label = 'SNP_23']
      av [pos = '-5,-5!', label = 'SNP_24']
     

      w->x
      a->w
      b->w
      c->w
      d->w
      e->w
      f->w
      g->w
      h->w
      i->w
      j->w
      k->w
      l->w
      m->w
      n->w
      o->w
      p->w
      q->w
      r->w
      s->w
      t->w
      u->w
      v->w
      y->x
      z->x
      az->x
      vv->w
   
      aw->x
      a1->aw
      aa->aw 
      ab->aw 
      ac->aw 
      ad->aw
      ae->aw
      af->aw
      ag->aw
      ah->aw
      ai->aw 
      aj->aw 
      ak->aw 
      al->aw 
      am->aw
      an->aw
      ao->aw
      ap->aw
      aq->aw
      ar->aw 
      as->aw
      at->aw
      au->aw
      av->aw
      a3->aw
      }
      ")


my_graph %>%export_svg %>% charToRaw %>% rsvg_png("Fig4_path_diag.png")