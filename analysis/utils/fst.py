import numpy as np

def calc_simple_fst_power(pops_gt, pops): 
    
    nsnps = pops_gt[pops[0]].shape[0]
    for j in range(1, len(pops)):
        if nsnps != pops_gt[pops[j]].shape[0]:
            print("SNP numbers differ between populations")
            raise
    
    power_list = list()
    fst_list = list()
    # all_gt = np.hstack((eur_gt1kg, afr_gt1kg))
    all_gt = np.hstack(tuple([pops_gt[x] for x in pops]))
    all_n  = all_gt.shape[1] # size of total population
    for snpi in range(all_gt.shape[0]):
        # maf_all = sum(all_gt[snpi,:] / 2 / len(all_gt[snpi,:]))
        # maf_eur = sum(eur_gt[snpi,:] / 2 / len(eur_gt[snpi,:]))
        # maf_afr = sum(afr_gt[snpi,:] / 2 / len(afr_gt[snpi,:]))
        if np.all(all_gt[snpi,:] == all_gt[snpi,0]):
            print("Warning!")
            Fst = -1.
            power = -1.
        else:
            maf_all = sum(all_gt[snpi,:] / 2 / len(all_gt[snpi,:]))
            maf_pops = [sum(pops_gt[x][snpi,:] / 2 / len(pops_gt[x][snpi,:])) for x in pops]
            
            # maf_all = sum(maf_pops) / len(maf_pops)

            c_pops = [pops_gt[x].shape[1]/all_n for x in pops]
            #c_eur = eur_n / all_n
            #c_afr = afr_n / all_n
            
            p_all = maf_all*(1-maf_all)
            sum_elems = np.array([c_pops[i]*maf_pops[i]*(1-maf_pops[i]) for i in range(len(pops))])
            sum_p = np.sum(sum_elems)
            power = all_n * sum_p # (c_eur*maf_eur*(1-maf_eur) + c_afr*maf_afr*(1-maf_afr) )
            
            Fst = (p_all  - sum_p ) / p_all
        power_list.append(power)
        fst_list.append(Fst)
    return fst_list, power_list



def get_gt_count(dosage, allele='ref'):
    n_het = dosage.count(1)
    if allele == "ref":
        n_hom = dosage.count(0)
    if allele == "alt":
        n_hom = dosage.count(2)
    return n_hom, n_het

def weir_cockerman_fst(pops_gt, pops=["eur", "afr"], alleles=["ref", "alt"]):
    
    fst_list = list()
    nsnps = pops_gt[pops[0]].shape[0]
    if nsnps != pops_gt[pops[1]].shape[0]:
        print("SNP numbers differ between populations")
        raise
        
    for snpi in range(nsnps):
        n_pops = len(pops)
        n_alleles = len(alleles)
        n = np.zeros(n_pops)
        p = np.zeros((n_pops, n_alleles))
        pbar = np.zeros(n_alleles)
        hbar = np.zeros(n_alleles)
        ssqr = np.zeros(n_alleles)

        nbar = 0.0
        sum_nsqr = 0.0
        n_sum = 0.0
        for pop in range(n_pops):
            for al in range(n_alleles):
                n_hom, n_het = get_gt_count(list(pops_gt[pops[pop]][snpi,:]), allele=alleles[al])
                n[pop] += n_hom + 0.5*n_het
                p[pop,al] = n_het + 2*n_hom

                nbar += n[pop]
                pbar[al] += p[pop][al]
                hbar[al] += n_het
            for al in range(n_alleles):
                p[pop,al] /= 2.0*n[pop]

            sum_nsqr += (n[pop] * n[pop])

        n_sum = sum(n)
        nbar  = n_sum / n_pops

        for al in range(n_alleles):
            pbar[al] /= (n_sum * 2.0)
            hbar[al] /= n_sum

        for al in range(n_alleles):
            for pop in range(n_pops):
                ssqr[al] += n[pop]*(p[pop,al] - pbar[al])*(p[pop,al] - pbar[al])
            ssqr[al] /= (n_pops-1)*nbar
        nc = (n_sum - (sum_nsqr / n_sum)) / (n_pops - 1)

        snp_Fst = np.zeros(n_alleles)
        a = np.zeros(n_alleles)
        b = np.zeros(n_alleles)
        c = np.zeros(n_alleles)
        r = n_pops
        sum_a = 0.0
        sum_all = 0.0
        for al in range(n_alleles):
            a[al] = (ssqr[al] - ( pbar[al]*(1.0-pbar[al]) - (((r-1.0)*ssqr[al])/r) - (hbar[al]/4.0) )/(nbar-1.0))*nbar/nc;
            b[al] = (pbar[al]*(1.0-pbar[al]) - (ssqr[al]*(r-1.0)/r) - hbar[al]*( ((2.0*nbar)-1.0) / (4.0*nbar) ))*nbar / (nbar-1.0) ;
            c[al] = hbar[al] / 2.0;
            snp_Fst[al] = a[al]/(a[al]+b[al]+c[al]);

            if not np.all([np.isnan(a[al]),np.isnan(b[al]),np.isnan(c[al])]):
                sum_a += a[al]
                sum_all += a[al] + b[al] + c[al]
        fst = sum_a/sum_all
        fst_list.append(fst)
    return fst_list