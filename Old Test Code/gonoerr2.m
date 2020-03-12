/cm/local/apps/torque/4.2.4.1/spool/mom_priv/jobs/9203648.master.cm.cluster.SC: line 1: fg: no job control
/cm/local/apps/torque/4.2.4.1/spool/mom_priv/jobs/9203648.master.cm.cluster.SC: line 2: fg: no job control
/cm/local/apps/torque/4.2.4.1/spool/mom_priv/jobs/9203648.master.cm.cluster.SC: line 9: PBS: command not found
{Undefined variable "VacAMR_IBM2" or class "VacAMR_IBM2.plfit".

Error in VacAMR_IBM3.pl_network (line 2169)
                [alpha_cont,x_cont,L_cont] = VacAMR_IBM2.plfit(node_degree);

Error in VacAMR_IBM3 (line 463)
                            [self.adj_full,~,self.rel_list_full] =
                            self.pl_network(self.N,self.ALPHA,self.FULL_MAX_PARTNERS);

Error in GonoTest (line 40)
        gono_model = VacAMR_IBM3(N, params, [], VERBOSE, LOW_MEM, vac); %turn
        on and off strategies here
} 
Warning: Permanently added 'node32-011,10.131.0.155' (RSA) to the list of known hosts.
