<!-- Sohan Das , CrS2119 , ISI,Kolkata -->

<?php
    $host="localhost:3307";
    $user="root";
    $password="";
    $db="sql_injection";

    $con=new mysqli($host,$user,$password,$db);

    if(!$con)
        echo "Could not connect with <b> Server </b>";
    else
        echo "Connected with <b> Server </b>";

    // mysql_select_db($db);

    if(isset($_GET['username'])){
            $uname=$_GET['username'];
            $password=$_GET['password'];

            $qry="select * from log_in where User='$uname' AND Password='$password'";

            $result=mysqli_query($con,$qry);

            

            if(mysqli_num_rows($result)){  
                echo "<br>";
                echo nl2br("<center>\n\n You Have Successfully <u> <b>Logged In</b></u></center>");

                header('Location: index.html');

                exit();
            }
            else{
                echo nl2br("<center>\n\n You Have Entered <u> <b>Wrong Inputs</b> </u></center>");
                exit();
            }
        }
?>



<!DOCTYPE html>
<html lang="en">

<head>
    <meta charset="UTF-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SQL_Injection_Attack</title>

    <script>
            function GEEKFORGEEKS() {
                var name =
                    document.forms.logForm.username.value;
                
                var regName = /^[
                    -zA-Z]+$/;                                    // Javascript reGex for Name validation
 
                if (name == "" || !regName.test(name)) {
                    window.alert("Please enter your name properly.");
                    name.focus();
                    return false;
                }
                return true;
            }

            const form = document.querySelector('#logForm');
            form.addEventListener('submit', function (e) {
                // prevent the form from submitting
                e.preventDefault();
                
            });  

    </script>


</head>

<body>
    <div style="text-align:center;color:red;">
        <h1> Login Form </h1>
    </div>
    
    <center>
    <div class="container" 
    style="width: 300px;
        height: 100px;
        padding: 70px;
        border: 3px solid black;
        background-color:LightSkyBlue;">
         <form id= logForm name=logForm method="GET" action="" style="text-align: center" onsubmit="return GEEKFORGEEKS()"> 
         <!-- onsubmit="return GEEKFORGEEKS()" -->
              <div class="form_input">
                  Username : <input type="text" id="username" name="username" placeholder="Enter Your Username" >
              </div>
              <br>
              <div class="form_input">
                  Password : <input type="text" name="password" placeholder="Enter Your Password">
              </div>
              <br>
              <input type="submit" name="submit" value="LOGIN"/>
              <br>
              <br>
              <a href="upPass.php">Forget password</a>
         </form>
    </div>
    </center>
</body>
</html>
    
     

    