
/**
 * Write MCMC results to parameter and model files.
 */
int io_write_results(IOConfig *ioconfig, MCMCResults *results) {
    FILE *fp;
    int nloci = ioconfig->model.marker_set_size; /* Wait, need nloci from somewhere. IOConfig doesn't track it directly? */
    /* Using genomic data size? No, I don't have gdata here. */
    /* I should assume results arrays are sized correctly. But I need to know the size loop. */
    /* Actually IOConfig doesn't store nloci. I need to pass it or rely on loading it. */
    /* Let's double check if I can get nloci from IOConfig. No. */
    /* I should probably pass nloci to this function or update IOConfig to store it after loading. */
    
    /* Better approach: Pass GenomicData or just count from file? */
    /* The caller (main) knows nloci. I should pass it. */
    /* But the signature I defined is (IOConfig*, MCMCResults*). MCMCResults doesn't store nloci size. */
    /* MCMCResults contains pointers. */
    /* I should update MCMCResults to include dimensions (nloci, nind) or pass them. */
    /* MCMCState/GenomicData had them. */
    /* Let's look at MCMCResults again. Use view_file. */

    return SUCCESS; 
}
