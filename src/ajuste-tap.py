import numpy as np
from cmath import polar, atan
from network import Network
from quadripoles import Quadripole, SeriesImpedance, ShuntAdmittance, PiCircuit, Transformer

vcondical3 = 69E3
vcondical2 = 230E3
vcondical1 = 500E3

aux = []
aux1 = []
aux2 = []

# Considera-se uma variável por vez

def findn1n2n3():
        # n1=3.33
        # n2=0.1
        n3=0.257

        while(1):

            network = Network(v1=69e3)
            s_imped_th = SeriesImpedance(z=4+0.38j)
            t1 = Transformer(n= 8.02, z1= 7.6e-3 + 3.8e-3j, z2= 33.9e-3 + 0.85e-3j, y=(4320 + 5050j)/(4320*5050j))

            network.set_resultant_matrix(network.cascade_connection(s_imped_th.net_mtrx, t1.net_mtrx))

            rkm = 0.172
            xlkm = 120*np.pi*2.18e-3j
            xckm = 1/(120*np.pi*0.0136e-6j)

            lt1 = PiCircuit(z= rkm*100 + xlkm*100, y1=1/(2*xckm) * 100, y2=1/(2*xckm) * 100)
            lt2 = PiCircuit(z= rkm*100 + xlkm*100, y1=1/(2*xckm) * 100, y2=1/(2*xckm) * 100)
            lt3 = PiCircuit(z= rkm*100 + xlkm*100, y1=1/(2*xckm) * 100, y2=1/(2*xckm) * 100)

            ltpl1 = network.parallel_connection(network.parallel_connection(lt1.net_mtrx, lt2.net_mtrx), lt3.net_mtrx)

            network.set_resultant_matrix(network.cascade_connection(network.resultant_mtrx, ltpl1))

            z1=ShuntAdmittance(z=8530 + 120*np.pi*52j)

            network.set_resultant_matrix(network.cascade_connection(network.resultant_mtrx, z1.net_mtrx))

            lt4 = PiCircuit(z= rkm*100 + xlkm*100, y1=1/(2*xckm) * 100, y2=1/(2*xckm) * 100)
            lt5 = PiCircuit(z= rkm*100 + xlkm*100, y1=1/(2*xckm) * 100, y2=1/(2*xckm) * 100)

            lt2pl = network.parallel_connection(lt4.net_mtrx, lt5.net_mtrx)

            network.set_resultant_matrix(network.cascade_connection(network.resultant_mtrx, lt2pl))

            t2 = Transformer(n= 0.5, z1= 7.6e-3 + 3.8e-3j, z2= 33.9e-3 + 0.85e-3j,
                            y= (432000 + 505000j)/(432000*505000j))

            network.set_resultant_matrix(network.cascade_connection(network.resultant_mtrx, t2.net_mtrx))

            z2=ShuntAdmittance(z=1050.55 + 120*np.pi*6.02j)


            network.set_resultant_matrix(network.cascade_connection(network.resultant_mtrx, z2.net_mtrx))

            lt6=PiCircuit(z= rkm*120 + xlkm*120, y1=1/(2*xckm) * 120, y2=1/(2*xckm) * 120)

            network.set_resultant_matrix(network.cascade_connection(network.resultant_mtrx, lt6.net_mtrx))

            t3=Transformer(n=n3, z1= 7.6e-3 + 3.8e-3j, z2= 33.9e-3 + 0.85e-3j,
                            y=(402000 + 607000j)/(402000*607000j))

            network.set_resultant_matrix(network.cascade_connection(network.resultant_mtrx, t3.net_mtrx))

            z3=ShuntAdmittance(z=505 + 120*np.pi*2.7j)


            def fasor(nmcomplexo):
                modulo, fase = polar(nmcomplexo)
                fase = np.rad2deg(fase)
                return modulo, fase

            vz3 = 69e3/(network.resultant_mtrx[0][0] + network.resultant_mtrx[0][1]/(505 + 120*np.pi*2.7j))
    
            iz3 = vz3/(505 + 120*np.pi*2.7j)

            fasorvz3 = fasor(vz3)
            fasoriz3 = fasor(iz3)

            net2 = Network()

            net2.set_resultant_matrix(net2.cascade_connection(z2.net_mtrx, lt6.net_mtrx))
            net2.set_resultant_matrix(net2.cascade_connection(net2.resultant_mtrx, t3.net_mtrx))

            vz2 = net2.resultant_mtrx[0][0]*vz3 + net2.resultant_mtrx[0][1]*iz3
            IS2  = net2.resultant_mtrx[1][0]*vz3 + net2.resultant_mtrx[1][1]*iz3

            iz2 = vz2/(1050.55 + 120*np.pi*6.02j)

            fasorvz2 = fasor(vz2)
            fasoriz2 = fasor(iz2)

            net3 = Network()

            net3.set_resultant_matrix(net3.cascade_connection(lt2pl, t2.net_mtrx))

            vz1 = net3.resultant_mtrx[0][0]*vz2 + net3.resultant_mtrx[0][1]*IS2
            iz1 = vz1/(8530 + 120*np.pi*52j)

            fasorvz1 = fasor(vz1)
            fasoriz1 = fasor(iz1)

            aux.append(fasorvz1[0])
            aux1.append(fasorvz2[0])
            aux2.append(fasorvz3[0])
            
            with open('dadosn1n2n3.txt','at') as f:
                f.write(f'n1:{8.02}, Tensao1: {fasorvz1[0]}, {fasorvz1[1]}\n')
                f.write(f'n2:{0.5}, Tensao2: {fasorvz2[0]}, {fasorvz2[1]}\n')
                f.write(f'n3:{n3}, Tensao3: {fasorvz3[0]}, {fasorvz3[1]}\n')
                f.write('-='*32)
                f.write('\n')
                
            # if fasorvz1[0] < 500e3:
            #     n1 += 0.1
            # if fasorvz2[0] < 230e3:
            #     n2 += 0.1
            if fasorvz3[0] < 69e3:
                n3 += 0.01
    
findn1n2n3()
