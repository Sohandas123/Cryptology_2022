Just run the program. All things are set. Also if you wish you can change the secret values  a  of Alice's or  b  of Bob's. We just check whether a secret private key is generated or not. If   a*(b*P) == b*(a*P)  then we are succeed and secure private key  abP is generated and they can use it as for further coomunications.  

Run the following command :
1.  gcc Elliptic_Curve_Key_Exchange.c functions.c -o output
2.  ./output