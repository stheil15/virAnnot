-- MySQL dump 10.13  Distrib 5.5.31, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: taxonomy
-- ------------------------------------------------------
-- Server version	5.5.31-0ubuntu0.12.04.2

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `gi_taxid_nucl`
--
PRAGMA encoding = "UTF-8";
PRAGMA default_synchronous = OFF;
CREATE TABLE `gi_taxid_nucl` (
  `gi` double NOT NULL,
  `tax_id` double NOT NULL,
  PRIMARY KEY (`gi`)
) ;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gi_taxid_prot`
--


CREATE TABLE `gi_taxid_prot` (
  `gi` double NOT NULL,
  `tax_id` double NOT NULL,
  PRIMARY KEY (`gi`)
) ;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `names`
--


CREATE TABLE `names` (
  `tax_id` int unsigned NOT NULL DEFAULT '0',
  `name_txt` varchar(255) NOT NULL DEFAULT '',
  `unique_name` varchar(255) DEFAULT NULL,
  `name_class` varchar(32) NOT NULL DEFAULT ''
) ;
/*!40101 SET character_set_client = @saved_cs_client */;


--
-- Table structure for table `nodes`
--


CREATE TABLE `nodes` (
  `tax_id` int unsigned NOT NULL DEFAULT '0',
  `parent_tax_id` int unsigned NOT NULL DEFAULT '0',
  `rank` varchar(32) DEFAULT NULL,
  `embl_code` varchar(16) DEFAULT NULL,
  `division_id` int NOT NULL DEFAULT '0',
  `inherited_div_flag` int NOT NULL DEFAULT '0',
  `genetic_code_id` int NOT NULL DEFAULT '0',
  `inherited_GC_flag` int NOT NULL DEFAULT '0',
  `mitochondrial_genetic_code_id` int NOT NULL DEFAULT '0',
  `inherited_MGC_flag` int NOT NULL DEFAULT '0',
  `GenBank_hidden_flag` int NOT NULL DEFAULT '0',
  `hidden_subtree_root_flag` int NOT NULL DEFAULT '0',
  `comments` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`tax_id`)
) ;


/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-07-02 10:17:11
