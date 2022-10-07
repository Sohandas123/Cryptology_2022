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

    if(isset($_POST['username'])){
            $uname=$_POST['username'];
            $password=$_POST['old_password'];
            $new_pass=$_POST['new_password'];

            $qry="update log_in set Password='$new_pass' 
            where User='$uname' AND Password='$password'";
            $result=mysqli_query($con,$qry);

            $chqry="select * from log_in where User='$uname' AND Password='$new_pass'";
            $rslt=mysqli_query($con,$chqry);

            if(mysqli_num_rows($rslt)){  
                echo "<br>";
                echo nl2br("<center>\n\n You Have Successfully <u> <b>UPDATED </b></u> Your Password </center>");

                // header('Location: index.html');

                exit();
            }
            else{
                echo nl2br("<center>\n\n You Have Entered <u> <b>Wrong Inputs</b> </u> Password is not Updated.
                </center>");
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
    <title>SQL_Injection_Update_Password</title>

</head>

<body>
    <div style="text-align:center;color:red;">
        <h1> Update Password Form </h1>
    </div>
    
    <center>
    <div class="container" 
    style="width: 300px;
        height: 100px;
        padding: 70px;
        border: 3px solid black;
        background-color:LightSkyBlue;">
         <form method="POST" action="" style="text-align: center">
              <div class="form_input">
              Username <input type="text" name="username" placeholder="Enter Your Username" >
              </div>
              <br>
              <div class="form_input">
              Old Password <input type="text" name="old_password" placeholder="Enter Your Old Password">
              </div>
              <br>
              <div class="form_input">
              New Password <input type="text" name="new_password" placeholder="Enter Your New Password">
              </div>
              <br> 
              <a href="login.php" class="Back">&laquo; Back to LOGIN</a> &nbsp &nbsp
              <input type="submit" name="submit" value="UPDATE">
              <br>

         </form>
    </div>
    </center>
</body>
</html>