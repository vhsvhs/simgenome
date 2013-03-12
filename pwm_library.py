"""This library contains functions that return DNA-binding specificities for several
transcription regulator genes.  The specifities are returned as multi-line strings,
with specifity values ranging from 0.0 to 1.0, and summing to 1.0 for each site in each motif."""


def get_random_pwm(n):
    """Returns a random PWM with n sites."""
    line = ""
    for i in range(0, n):
        vals = []
        sum = 0.0
        for i in range(0,4):
            x = random.random()
            vals.append( x )
            sum += x
        for val in vals:
            line += (val/sum).__str__() + " "
        line += "\n"
    return line

def get_nonspec_pwm(n):
    """Returns a totally non-specific (all specificities = 0.25), with n sites."""
    line = ""
    for i in range(0, n):
        line += "0.25 0.25 0.25 0.25\n"
    return line

def get_Aspec(n):
    """Returns an A-specific PWM with n sites."""
    line = ""
    for i in range(0, n):
        line += "1.00 0.0 0.0 0.0\n"
    return line

def get_Cspec(n):
    """Returns an C-specific PWM with n sites."""
    line = ""
    for i in range(0, n):
        line += "0.00 1.0 0.0 0.0\n"
    return line

def get_Gspec(n):
    """Returns an G-specific PWM with n sites."""
    line = ""
    for i in range(0, n):
        line += "0.00 0.0 1.0 0.0\n"
    return line

def get_Tspec(n):
    """Returns an T-specific PWM with n sites."""
    line = ""
    for i in range(0, n):
        line += "0.00 0.0 0.0 1.0\n"
    return line

def get_matalpha2():
    # ATTTACATG
    line = ""
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.60 0.01 0.01 0.38\n"
    line += "0.10 0.01 0.09 0.80\n"
    line += "0.01 0.01 0.01 0.97\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.01 0.97 0.01 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.01 0.20 0.01 0.78\n"
    line += "0.08 0.01 0.90 0.01\n"    
    return line

def get_ste12():
    line = ""
    line += "0.01 0.01 0.01 0.97\n"
    line += "0.01 0.01 0.97 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.01 0.97 0.01 0.01\n"
    line += "0.70 0.01 0.28 0.01\n"  
    return line

def get_tec1():
    """Returns the TEC1 PWM, from Nobile et al., 2012"""
    line = ""
    line += "0.201541    0.113628    0.562196    0.122633\n"
    line += "0.377366    0.014727    0.606152    0.001754\n"
    line += "0.025716    0.036705    0.935823    0.001754\n"
    line += "0.981764    0.014727    0.001754    0.001754\n"
    line += "0.992753    0.003738    0.001754    0.001754\n"
    line += "0.014727    0.014727    0.001754    0.968791\n"
    line += "0.003738    0.003738    0.97978     0.012743\n"
    line += "0.080661    0.223519    0.056699    0.639119\n"
    line += "0.377366    0.102639    0.100655    0.419338\n"
    return line

def get_bcr1():
    """Returns the BCR1 PWM, from Nobile et al., 2012."""
    line = ""
    line += "0.026635    0.026635    0.011805    0.934924\n"
    line += "0.949754    0.026635    0.011805    0.011805\n"
    line += "0.180488    0.795901    0.011805    0.011805\n"
    line += "0.949754    0.026635    0.011805    0.011805\n"
    line += "0.026635    0.026635    0.011805    0.934924\n"
    line += "0.334341    0.026635    0.627217    0.011805\n"
    line += "0.026635    0.949754    0.011805    0.011805\n"
    line += "0.949754    0.026635    0.011805    0.011805\n"
    line += "0.026635    0.180488    0.011805    0.78107\n"
    line += "0.565121    0.103561    0.319511    0.011805\n"
    line += "0.411268    0.257415    0.011805    0.319511\n"
    line += "0.642047    0.257415    0.011805    0.088731\n"
    return line 

def get_brg1():
    """Returns the BRG1 PWM, from Nobile et al., 2012."""
    line = ""
    line += "0.031478    0.486048    0.468521    0.013951\n"
    line += "0.75879        0.213306    0.013951    0.013951\n"
    line += "0.031478    0.031478    0.923091    0.013951\n"
    line += "0.031478    0.031478    0.923091    0.013951\n"
    line += "0.031478    0.031478    0.013951    0.923091\n"
    line += "0.940618    0.031478    0.013951    0.013951\n"
    line += "0.213306    0.75879        0.013951    0.013951\n"
    return line

def get_efg1():
    """Returns the EFG1 PWM, from Nobile et al., 2012."""
    line = ""
    line += "0.164769    0.140078    0.026586    0.668566\n"
    line += "0.01662        0.547488    0.001894    0.433996\n"
    line += "0.954898    0.01662        0.026586    0.001894\n"
    line += "0.01662        0.01662        0.001894    0.964864\n"
    line += "0.033241    0.008549    0.929728    0.02848\n"
    line += "0.004274    0.991936    0.001894    0.001894\n"
    line += "0.868478    0.01662        0.10066        0.01424\n"
    line += "0.041311    0.522796    0.01424        0.42165\n"
    return line

def get_ndt80():
    """Returns the NDT80 PWM, from Nobile et al., 2012"""
    line = ""
    line += "0.016417    0.028612    0.123823    0.831145\n"
    line += "0.028612    0.016417    0.245775    0.709193\n"
    line += "0.882277    0.016417    0.038457    0.062847\n"
    line += "0.016417    0.967644    0.001871    0.014066\n"
    line += "0.979839    0.016417    0.001871    0.001871\n"
    line += "0.016417    0.979839    0.001871    0.001871\n"
    line += "0.992034    0.004222    0.001871    0.001871\n"
    line += "0.992034    0.004222    0.001871    0.001871\n"
    line += "0.955448    0.016417    0.001871    0.026261\n"
    line += "0.857887    0.065198    0.050652    0.026261\n"
    line += "0.272517    0.540811    0.038457    0.148214\n"
    line += "0.199345    0.577397    0.038457    0.184799\n"   
    return line

def get_rob1():
    """Returns the ROB1 PWM, from Nobile et al., 2012."""
    line = ""
    line += "0.036962    0.036962    0.907489    0.018585\n"
    line += "0.036962    0.036962    0.907489    0.018585\n"
    line += "0.814753    0.036962    0.018585    0.129698\n"
    line += "0.925866    0.036962    0.018585    0.018585\n"
    line += "0.814753    0.148075    0.018585    0.018585\n"
    line += "0.259188    0.148075    0.018585    0.57415\n"
    line += "0.259188    0.148075    0.240811    0.351924\n"
    line += "0.70364        0.148075    0.018585    0.129698\n"
    line += "0.925866    0.036962    0.018585    0.018585\n"
    line += "0.370301    0.148075    0.018585    0.463037\n"
    line += "0.259188    0.036962    0.129698    0.57415\n"
    line += "0.036962    0.036962    0.018585    0.907489\n"
    line += "0.036962    0.925866    0.018585    0.018585\n"
    line += "0.036962    0.925866    0.018585    0.018585\n"
    line += "0.370301    0.036962    0.57415        0.018585\n"
    line += "0.370301    0.592527    0.018585    0.018585\n"
    return line
