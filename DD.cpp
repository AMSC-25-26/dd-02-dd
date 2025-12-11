/* LA CLASSE CHE CREA L'OGGETTO PER RISOLVERE IL PROBLEMA */  

class SchwarzSubdomain {

public:
    SchwarzSubdomain(int rank, int n_procs, double mu_val, double c_val, 
                     unsigned int n_elems, double overlap)
        : mpi_rank(rank), mpi_size(n_procs), 
          n_elements_local(n_elems), delta(overlap),
          mu(mu_val), c(c_val), fe(1),                
          dof_handler(triangulation) 
    {}

    void run(); 

private:

        /* 
        QUESTI SONO I METODI PER LA RISOLUZIONE DEL PROBLEMA, LI METTO PRIVATI PERCHE DEVO UTILIZZARLI SOLO IN run()
        1. CREO LA GRIGLIA
        2. INIZIALIZZA LE MATRICI E I VETTORI PER PREPARARE IL METODO ITERATIVO
        3. IMPONGO LE BOUNDARY CONDITIONS
        4. IMPLEMENTO IL METODO ITERATIVO (CG)
        5. METODO PER AVERE I RISULTATI IN FORMATO .vtU
        6. PRENDO I DATI CALCOLATI
        7. IMPONGO UN MODO PER FAR "COMUNICARE" DUE DOMINI AL BORDO
        */

    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results(int iteration);
    
    double get_solution_at_point(double x); 
    void communicate_boundaries();


    /* OGGETTI DI TEMPLATE DI deal.ii UTILI 
    PER CREARE IL PROBLEMA E RISOLVERLO CON GLI ELEMENTI FINITI*/

    Triangulation<1>     triangulation;
    FE_Q<1>              fe; 
    DoFHandler<1>        dof_handler;
    


    AffineConstraints<double> constraints; 
    
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs;
    Vector<double>       solution;
};
